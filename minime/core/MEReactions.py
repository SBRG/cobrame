from __future__ import division

from warnings import warn
from collections import defaultdict

from six import iteritems

from cobra import Reaction

from minime.util.mass import *
from minime.util import mu, dogma
from minime.core.Components import *

from ecolime.ribosome import translation_stop_dict


class MEReaction(Reaction):
    def add_modifications(self, process_data_id, stoichiometry):
        all_modifications = self._model.modification_data
        process_info = self._model.process_data.get_by_id(process_data_id)
        for modification_id, count in iteritems(process_info.modifications):
            modification = all_modifications.get_by_id(modification_id)
            for mod_comp, mod_count in iteritems(modification.stoichiometry):
                stoichiometry[mod_comp] += count * mod_count

            if type(modification.enzyme) == list:
                for enzyme in modification.enzyme:
                    stoichiometry[enzyme] -= \
                        mu / modification.keff / 3600.
            elif type(modification.enzyme) == str:
                stoichiometry[modification.enzyme] -= \
                    mu / modification.keff / 3600.

        return stoichiometry

    def add_subreactions(self, process_data_id, stoichiometry):
        all_subreactions = self._model.subreaction_data
        process_info = self._model.process_data.get_by_id(process_data_id)
        for subreaction_id, count in iteritems(process_info.subreactions):
            subreaction_data = all_subreactions.get_by_id(subreaction_id)

            if type(subreaction_data.enzyme) == list:
                for enzyme in subreaction_data.enzyme:
                    stoichiometry[enzyme] -= mu / subreaction_data.keff / 3600.
            elif type(subreaction_data.enzyme) == str:
                stoichiometry[subreaction_data.enzyme] -= \
                    mu / subreaction_data.keff / 3600.

            for met, stoich in subreaction_data.stoichiometry.items():
                stoichiometry[met] += count * stoich

        return stoichiometry

    def get_components_from_ids(self, id_stoichiometry,
                                default_type=Component, verbose=True):
        stoic = id_stoichiometry
        object_stoichiometry = {}
        mets = self._model.metabolites
        for key, value in iteritems(stoic):
            try:
                object_stoichiometry[mets.get_by_id(key)] = value
            except KeyError:
                new_met = create_component(key, default_type=default_type)
                if verbose:
                    print("Created %s in %s" %
                          (repr(new_met), repr(self)))
                object_stoichiometry[new_met] = value
                self._model.add_metabolites([new_met])
        return object_stoichiometry


class MetabolicReaction(Reaction):
    """Metabolic reaction including required enzymatic complex"""

    @property
    def complex_data(self):
        return self._complex_data

    @complex_data.setter
    def complex_data(self, process_data):
        self._complex_data = process_data
        process_data._parent_reactions.add(self.id)
    _complex_data = None

    @property
    def stoichiometric_data(self):
        return self._stoichiometric_data

    @stoichiometric_data.setter
    def stoichiometric_data(self, process_data):
        self._stoichiometric_data = process_data
        process_data._parent_reactions.add(self.id)
    _stoichiometric_data = None

    keff = 65.  # in per second
    reverse = False
    complex_dilution_list = set()

    def update(self, create_new=False):
        metabolites = self._model.metabolites
        new_stoichiometry = defaultdict(float)
        if self.complex_data:
            new_stoichiometry[self.complex_data.complex] = \
                -mu / self.keff / 3600.
        sign = -1 if self.reverse else 1
        for component_id, value in iteritems(
                self.stoichiometric_data._stoichiometry):
            if create_new:
                try:
                    component = metabolites.get_by_id(component_id)
                except KeyError:
                    component = create_component(component_id)
                    self._model.add_metabolites([component])
                    print("Created %s in %s" % (repr(component), repr(self)))
            else:
                component = metabolites.get_by_id(component_id)
            if component not in new_stoichiometry:
                new_stoichiometry[component] = 0
            new_stoichiometry[component] += value * sign
            if component.id in self.complex_dilution_list:
                new_stoichiometry[component] += - mu / self.keff / 3600.
        # replace old stoichiometry with new one
        # TODO prune out old metabolites
        # doing checks for relationship every time is unnecessary. don't do.
        self.add_metabolites(new_stoichiometry,
                             combine=False, add_to_container_model=False)

        # set the bounds
        if self.reverse:
            self.lower_bound = max(
                0, -self.stoichiometric_data.upper_bound)
            self.upper_bound = max(
                0, -self.stoichiometric_data.lower_bound)
        else:
            self.lower_bound = max(0, self.stoichiometric_data.lower_bound)
            self.upper_bound = max(0, self.stoichiometric_data.upper_bound)


class ComplexFormation(MEReaction):
    """Formation of a protein complex"""
    _complex_id = None

    @property
    def complex(self):
        return self._model.metabolites.get_by_id(self._complex_id)

    def update(self, verbose=True):
        metabolites = self._model.metabolites
        complex_info = self._model.complex_data.get_by_id(self._complex_id)
        try:
            complex_met = metabolites.get_by_id(self._complex_id)
        except KeyError:
            complex_met = create_component(self._complex_id,
                                           default_type=Complex)
            self._model.add_metabolites([complex_met])
        stoichiometry = defaultdict(float)
        stoichiometry[complex_met.id] = 1

        # build the complex itself
        for component_id, value in iteritems(complex_info.stoichiometry):
            stoichiometry[component_id] -= value

        # add in the modifications
        stoichiometry = self.add_modifications(self._complex_id, stoichiometry)

        object_stoichiometry = self.get_components_from_ids(
                stoichiometry, default_type=Complex, verbose=verbose)

        self.add_metabolites(object_stoichiometry, combine=False,
                             add_to_container_model=True)


class TranscriptionReaction(Reaction):
    """Transcription of a TU to produced TranscribedGene"""
    _transcription_data = None

    @property
    def transcription_data(self):
        return self._transcription_data

    @transcription_data.setter
    def transcription_data(self, process_data):
        self._transcription_data = process_data
        process_data._parent_reactions.add(self.id)

    def update(self):
        protein_id = self.transcription_data.id
        new_stoichiometry = defaultdict(int)
        TU_length = len(self.transcription_data.nucleotide_sequence)
        metabolites = self._model.metabolites

        if self.transcription_data.using_RNAP:
            try:
                RNAP = self._model.metabolites.get_by_id("RNA_Polymerase")
            except KeyError:
                warn("RNA Polymerase not found")
            else:
                k_RNAP = (mu * 22.7 / (mu + 0.391))*3  # 3*k_ribo
                coupling = -TU_length * mu / k_RNAP / 3600
                new_stoichiometry[RNAP] = coupling

        for transcript_id in self.transcription_data.RNA_products:
            if transcript_id not in self._model.metabolites:
                transcript = TranscribedGene(transcript_id)
                self._model.add_metabolites([transcript])
            else:
                transcript = metabolites.get_by_id(transcript_id)
            new_stoichiometry[transcript] += 1

        # add in the modifications
        all_modifications = self._model.modification_data
        for modification_id, count in iteritems(
                self.transcription_data.modifications):
            modification = all_modifications.get_by_id(modification_id)
            for mod_comp, mod_count in iteritems(modification.stoichiometry):
                new_stoichiometry[metabolites.get_by_id(mod_comp)] += count * mod_count
            if modification.enzyme is not None:
                new_stoichiometry[metabolites.get_by_id(modification.enzyme)] -= \
                    mu / modification.keff / 3600.

        base_counts = self.transcription_data.nucleotide_count

        for base, count in iteritems(base_counts):
            new_stoichiometry[metabolites.get_by_id(base)] -= count
        for base, count in iteritems(self.transcription_data.excised_bases):
            new_stoichiometry[metabolites.get_by_id(base)] += count

        new_stoichiometry[metabolites.get_by_id("ppi_c")] += TU_length - 1
        # 5' had a triphosphate, but this is removed when excising
        # Is this universally true?
        if sum(self.transcription_data.excised_bases.values()) > 0:
            new_stoichiometry[metabolites.get_by_id("ppi_c")] += 1

        self.add_metabolites(new_stoichiometry, combine=False,
                             add_to_container_model=False)

        # add to biomass
        RNA_mass = self.transcription_data.mass  # kDa
        self.add_metabolites({self._model._biomass: RNA_mass},
                             combine=False, add_to_container_model=False)


class GenericFormationReaction(Reaction):
    None


class TranslationReaction(MEReaction):
    """Translation of a TranscribedGene to a TranslatedGene"""
    _translation_data = None

    @property
    def translation_data(self):
        return self._translation_data

    @translation_data.setter
    def translation_data(self, process_data):
        self._translation_data = process_data
        process_data._parent_reactions.add(self.id)

    def update(self, verbose=True):
        protein_id = self.translation_data.protein
        mRNA_id = self.translation_data.mRNA
        protein_length = len(self.translation_data.amino_acid_sequence)
        metabolites = self._model.metabolites
        new_stoichiometry = defaultdict(int)

        # If ribosome is not added to model, warn user but do not raise error
        if self.translation_data.using_ribosome:
            try:
                ribosome = self._model.metabolites.get_by_id("ribosome")
            except KeyError:
                warn("ribosome not found")
            else:
                k_ribo = mu * 22.7 / (mu + 0.391)
                coupling = -protein_length * mu / k_ribo / 3600.
                new_stoichiometry[ribosome.id] = coupling

        # Not all genes have annotated TUs. For these add TU as the bnumber
        try:
            transcript = metabolites.get_by_id(mRNA_id)
        except KeyError:
            warn("transcript '%s' not found" % mRNA_id)
            transcript = TranscribedGene(mRNA_id)
            self._model.add_metabolites(transcript)

        k_mRNA = mu * self.translation_data.protein_per_mRNA / (mu + 0.391)
        new_stoichiometry[transcript.id] = -mu / k_mRNA / 3600.

        # Added protein to model if not already included
        try:
            protein = metabolites.get_by_id(protein_id)
        except KeyError:
            protein = TranslatedGene(protein_id)
            self._model.add_metabolites(protein)
        new_stoichiometry[protein.id] = 1

        # ------------------ Elongation Reactions------------------------
        # update stoichiometry
        aa_count = self.translation_data.amino_acid_count
        for aa_name, value in iteritems(aa_count):
            new_stoichiometry[aa_name] -= value

        # add in the tRNA's for each of the amino acids
        # We add in a "generic tRNA" for each amino acid. The production
        # of this generic tRNA represents the production of enough of any tRNA
        # (according to its keff) for that amino acid to catalyze addition
        # of a single amino acid to a growing peptide.

        # Add transcription termination as subprocessdata so keff can be
        # modified
        all_subreactions = self._model.subreaction_data
        if self.translation_data.using_ribosome:
            for codon, count in self.translation_data.codon_count.items():
                codon = codon.replace('U','T')
                if codon == 'TGA':
                    print 'Adding selenocystein for %s' % mRNA_id
                    aa = 'sec'
                else:
                    abbreviated_aa = dogma.codon_table[codon]
                    aa = dogma.amino_acids[abbreviated_aa].split('_')[0]
                codon = codon.replace('T', 'U')
                subreaction_id = aa + '_addition_at_' + codon
                try:
                    self.translation_data.subreactions[subreaction_id] = count
                except KeyError:
                    warn('subreaction %s not in model' % subreaction_id)

        self.translation_data.subreactions['fmet_addition_at_START'] = 1

        self.add_subreactions(self.translation_data.id, new_stoichiometry)


        # TODO: how many protons/water molecules are exchanged when making the
        # peptide bond?
        # tRNA charging requires 2 ATP per amino acid. TODO: tRNA demand
        # Translocation and elongation each use one GTP per event
        # (len(protein) - 1). Initiation and termination also each need
        # one GTP each. Therefore, the GTP hydrolysis resulting from
        # translation is 2 * len(protein)

        # tRNA + GTP -> tRNA_GTP
        new_stoichiometry["h2o_c"] -= 3 * protein_length
        new_stoichiometry["h_c"] += 3 * protein_length
        new_stoichiometry["pi_c"] += 3 * protein_length
        new_stoichiometry["gtp_c"] -= 2 * protein_length
        new_stoichiometry["gdp_c"] += 2 * protein_length
        new_stoichiometry["atp_c"] -= 1 * protein_length
        new_stoichiometry["adp_c"] += 1 * protein_length

        last_codon = \
            self.translation_data.nucleotide_sequence[-3:].replace('T', 'U')

        term_enzyme = translation_stop_dict.get(last_codon)
        try:
            term_subreaction_data = all_subreactions.get_by_id(
                last_codon + '_' + term_enzyme + '_mediated_termination')
            new_stoichiometry[term_subreaction_data.enzyme] -= \
                mu / term_subreaction_data.keff / 3600.
        except:
            warn('Term Enzyme not in model for %s' % self.id)

        object_stoichiometry = self.get_components_from_ids(new_stoichiometry,
                                                            verbose=verbose)

        self.add_metabolites(object_stoichiometry,
                             combine=False, add_to_container_model=False)

        # add to biomass
        protein_mass = compute_protein_mass(aa_count)  # kDa
        self.add_metabolites({self._model._biomass: protein_mass},
                             combine=False, add_to_container_model=False)


class tRNAChargingReaction(MEReaction):

    tRNAData = None

    def update(self, verbose=True):
        stoic = defaultdict(int)
        data = self.tRNAData
        mets = self._model.metabolites
        # If the generic tRNA does not exist, create it now. The meaning of
        # a generic tRNA is described in the TranslationReaction comments
        generic_tRNA = "generic_tRNA_" + data.codon + "_" + data.amino_acid

        stoic[generic_tRNA] = 1

        # compute what needs to go into production of a generic tRNA
        tRNA_amount = mu / data.tRNA_keff / 3600
        synthetase_amount = mu / data.synthetase_keff / \
                            3600 * (1 + tRNA_amount)
        stoic[data.RNA] = -tRNA_amount
        stoic[data.amino_acid] = -tRNA_amount
        if data.synthetase is not None:
            stoic[data.synthetase] = -synthetase_amount

        stoic = self.add_modifications(self.tRNAData.id, stoic)

        object_stoichiometry = self.get_components_from_ids(stoic,
                                                            verbose=verbose)

        self.add_metabolites(object_stoichiometry, combine=False,
                             add_to_container_model=False)


class SummaryVariable(Reaction):
    pass
