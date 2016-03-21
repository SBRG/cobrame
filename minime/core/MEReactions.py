from __future__ import division

from warnings import warn
from collections import defaultdict

from six import iteritems

from cobra import Reaction

from minime.util.mass import *
from minime.util import mu, dogma
from sympy import Basic
from minime.core.Components import *


class MEReaction(Reaction):

    # TODO set _upper and _lower bounds as a property
    """
    MEReaction is a general reaction class from which all ME-Model reactions
    will inherit

    This class contains functions used by all ME-model reactions


    """
    def add_modifications(self, process_data_id, stoichiometry, scale=1.):
        """
        Function to add modification process data to reaction stoichiometry


        process_data_id: String
            ID of the process data associated with the metabolic reaction.

            For example, if the modifications are being added to a complex
            formation reaction, the process data id would be the name of the
            complex.


        stoichiometry: Dictionary {metabolite_id: float} or
                                  {metabolite_id: float * (sympy.Symbol)}

        scale: float
           Some processes (ie. tRNA charging) are reformulated such other
           involved metabolites need scaling


        return: stoichiometry
            The dictionary with updated entries
        """

        all_modifications = self._model.modification_data
        process_info = self._model.process_data.get_by_id(process_data_id)
        for modification_id, count in iteritems(process_info.modifications):
            modification = all_modifications.get_by_id(modification_id)
            for mod_comp, mod_count in iteritems(modification.stoichiometry):
                stoichiometry[mod_comp] += count * mod_count * scale

            if type(modification.enzyme) == list:
                for enzyme in modification.enzyme:
                    stoichiometry[enzyme] -= \
                        mu / modification.keff / 3600.
            elif type(modification.enzyme) == str:
                stoichiometry[modification.enzyme] -= \
                    mu / modification.keff / 3600.

        return stoichiometry

    def add_subreactions(self, process_data_id, stoichiometry):
        """
        Function to add modification process data to reaction stoichiometry

        process_data_id: String
            ID of the process data associated with the metabolic reaction.

            For example, if the modifications are being added to a complex
            formation reaction, the process data id would be the name of the
            complex.


        stoichiometry: Dictionary {metabolite_id: float} or
                                  {metabolite_id: float * (sympy.Symbol)}


        return: stoichiometry
            The dictionary with updated entries
        """

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
        """
        Function to convert stoichiometry dictionary entries from strings to
        cobra objects.

        {metabolite_id: value} to {cobra.core.Metabolite: value}

        id_stoichiometry: Dict {string: float}
            Input Dict of {metabolite_id: value}

        default_type: String
            The type of cobra.Metabolite to default to if the metabolite is not
             yet present in the model

        verbose: Boolean
            If True, print metabolites added to model if not yet present in
            model

        return: object_stoichiometry
            {cobra.core.Metabolite: value}
        """

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

    def check_ME_mass_balance(self):
        """Compute mass and charge balance for the reaction

        returns a dict of {element: amount} for unbalanced elements.
        "charge" is treated as an element in this dict
        This should be empty for balanced reactions.
        """
        reaction_element_dict = defaultdict(int)
        for metabolite, coefficient in self._metabolites.items():
            if metabolite.id == 'biomass':
                continue
            if isinstance(coefficient, Basic):
                coefficient = coefficient.subs(mu, 0)
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += \
                    coefficient * metabolite.charge
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coefficient * amount
        # filter out 0 values
        return {k: v for k, v in iteritems(reaction_element_dict) if v != 0}


class MetabolicReaction(MEReaction):
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
    complex_data_id = None

    @property
    def complex(self):
        return self._model.metabolites.get_by_id(self._complex_id)

    def update(self, verbose=True):
        stoichiometry = defaultdict(float)
        metabolites = self._model.metabolites
        complex_info = self._model.complex_data.get_by_id(self.complex_data_id)
        try:
            complex_met = metabolites.get_by_id(self._complex_id)
        except KeyError:
            complex_met = create_component(self._complex_id,
                                           default_type=Complex)
            self._model.add_metabolites([complex_met])
        stoichiometry[complex_met.id] = 1

        # build the complex itself
        for component_id, value in iteritems(complex_info.stoichiometry):
            stoichiometry[component_id] -= value

        # add in the modifications
        stoichiometry = self.add_modifications(complex_info.id, stoichiometry)

        object_stoichiometry = self.get_components_from_ids(
                stoichiometry, default_type=Complex, verbose=verbose)

        # Add formula as sum of all protein and modification components
        elements = defaultdict(int)
        for component, value in iteritems(stoichiometry):
            if isinstance(value, Basic):
                value = value.subs(mu, 0)
            if component == complex_met.id:
                continue
            for e, n in iteritems(metabolites.get_by_id(component).elements):
                elements[e] += n * -int(value)

        complex_met.formula = "".join(stringify(e, n)
                                      for e, n in sorted(iteritems(elements)))

        self.add_metabolites(object_stoichiometry, combine=False,
                             add_to_container_model=True)


class TranscriptionReaction(MEReaction):
    """Transcription of a TU to produced TranscribedGene"""

    # TODO double check how initiation use is used as well as ATP cost etc.
    _transcription_data = None

    @property
    def transcription_data(self):
        return self._transcription_data

    @transcription_data.setter
    def transcription_data(self, process_data):
        self._transcription_data = process_data
        process_data._parent_reactions.add(self.id)

    def update(self, verbose=True):
        TU_id = self.transcription_data.id
        stoichiometry = defaultdict(int)
        TU_length = len(self.transcription_data.nucleotide_sequence)
        RNA_polymerase = self.transcription_data.RNA_polymerase
        metabolites = self.model.metabolites
        codes_stable_rna = self.transcription_data.codes_stable_rna  # TODO delete

        try:
            RNAP = self._model.metabolites.get_by_id(RNA_polymerase)
        except KeyError:
            warn("RNA Polymerase (%s) not found" % RNA_polymerase)
        else:
            k_RNAP = (mu * 22.7 / (mu + 0.391))*3  # 3*k_ribo
            coupling = -TU_length * mu / k_RNAP / 3600
            stoichiometry[RNAP.id] = coupling

        # This is necessary as opposed to using get_components_from_ids in
        # so that self.excised_bases can be used for "dummy_RNA"
        # TODO Can we fix this?
        for transcript_id in self.transcription_data.RNA_products:
            if transcript_id not in self._model.metabolites:
                transcript = TranscribedGene(transcript_id)
                self._model.add_metabolites([transcript])
            else:
                transcript = self._model.metabolites.get_by_id(transcript_id)
            stoichiometry[transcript.id] += 1

            # Add in formula for each transcript
            elements = defaultdict(int)
            for nucleotide, value in iteritems(transcript.nucleotide_count):
                for e, n in iteritems(metabolites.get_by_id(nucleotide).elements):
                    elements[e] += n * value

            transcript.formula = "".join(stringify(e, n)
                                         for e, n in sorted(iteritems(elements)))

        # add in the modifications
        stoichiometry = self.add_modifications(TU_id, stoichiometry)
        stoichiometry = self.add_subreactions(TU_id, stoichiometry)

        base_counts = self.transcription_data.nucleotide_count
        for base, count in iteritems(base_counts):
            stoichiometry[base] -= count
        for base, count in iteritems(self.transcription_data.excised_bases):
            stoichiometry[base] += count

        stoichiometry["ppi_c"] += TU_length - 1
        stoichiometry["h2o_c"] -= TU_length - 1
        stoichiometry["h_c"] += TU_length - 1
        # 5' had a triphosphate, but this is removed when excising
        # Is this universally true?
        if sum(self.transcription_data.excised_bases.values()) > 0:
            stoichiometry["ppi_c"] += 1

        new_stoich = \
            self.get_components_from_ids(stoichiometry, verbose=verbose,
                                         default_type=TranscribedGene)

        self.add_metabolites(new_stoich, combine=False,
                             add_to_container_model=False)

        # biomass contribution handled in translation
        RNA_mass = 0
        for met, v in new_stoich.items():
            if v < 0:
                continue
            if hasattr(met, "RNA_type") and met.RNA_type != 'mRNA':
                RNA_mass += met.formula_weight / 1000.  # kDa
        if RNA_mass > 0:
            self.add_metabolites({self._model._RNA_biomass: RNA_mass},
                                 combine=False, add_to_container_model=False)



class GenericFormationReaction(MEReaction):
    None


def stringify(element, number):
    return element if number == 1 else (element + str(number).rstrip("."))


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
        model = self._model
        metabolites = self._model.metabolites
        new_stoichiometry = defaultdict(int)

        try:
            ribosome = metabolites.get_by_id("ribosome")
        except KeyError:
            warn("ribosome not found")
        else:
            k_ribo = mu * 22.7 / (mu + 0.391)  # TODO move to translation_data
            coupling = -protein_length * mu / k_ribo / 3600.
            new_stoichiometry[ribosome.id] = coupling

        # Not all genes have annotated TUs. For these add TU as the bnumber
        try:
            transcript = metabolites.get_by_id(mRNA_id)
        except KeyError:
            warn("transcript '%s' not found" % mRNA_id)
            transcript = TranscribedGene(mRNA_id)
            model.add_metabolites(transcript)

        # Set Parameters
        kt = 4.5
        r0 = 0.087
        c_ribo = 1700./0.109/0.86
        m_aa = 0.109   # kDa
        k_mRNA = mu * self.translation_data.protein_per_mRNA / (mu + 0.391)
        RNA_amount = mu / k_mRNA / 3600.

        k_deg = 1.0/5. * 60.0  # 1/5 1/min 60 min/h # h-1 # TODO move to translation data
        deg_fraction = 3. * k_deg / (3. * k_deg + mu)
        deg_amount = deg_fraction * RNA_amount

        # set stoichiometry
        new_stoichiometry[transcript.id] = -RNA_amount
        for nucleotide, count in transcript.nucleotide_count.items():
            new_stoichiometry[nucleotide] += count * deg_amount
        nucleotide_length = len(transcript.nucleotide_sequence)
        # ATP hydrolysis required for cleaving
        new_stoichiometry['atp_c'] -= (nucleotide_length-1) / 4. * deg_amount
        new_stoichiometry['h2o_c'] -= (nucleotide_length-1) / 4. * deg_amount
        new_stoichiometry['adp_c'] += (nucleotide_length-1) / 4. * deg_amount
        new_stoichiometry['pi_c'] += (nucleotide_length-1) / 4. * deg_amount
        new_stoichiometry['h_c'] += (nucleotide_length-1) / 4. * deg_amount
        try:
            RNA_degradosome = metabolites.get_by_id("RNA_degradosome")
        except KeyError:
            warn("RNA_degradosome not found")
        else:
            new_stoichiometry[RNA_degradosome.id] = -deg_amount * mu / (65. * 3600.) #keff of degradosome

        # Added protein to model if not already included
        try:
            protein = metabolites.get_by_id(protein_id)
        except KeyError:
            protein = TranslatedGene(protein_id)
            model.add_metabolites(protein)
        new_stoichiometry[protein.id] = 1

        # ------------------ Elongation Reactions------------------------
        # update stoichiometry
        aa_count = self.translation_data.amino_acid_count
        elements = defaultdict(int)
        for aa_name, value in iteritems(aa_count):
            for e, n in iteritems(metabolites.get_by_id(aa_name).elements):
                elements[e] += n * value
        # subtract water from composition
        elements["H"] -= (protein_length - 1) * 2
        elements["O"] -= protein_length - 1
        # methionine becomes formylmethionine
        elements["C"] += 1
        elements["O"] += 1

        protein.formula = "".join(stringify(e, n)
                                  for e, n in sorted(iteritems(elements)))

        # add in the tRNA's for each of the amino acids
        # We add in a "generic tRNA" for each amino acid. The production
        # of this generic tRNA represents the production of enough of any tRNA
        # (according to its keff) for that amino acid to catalyze addition
        # of a single amino acid to a growing peptide.

        # Add transcription termination as subprocessdata so keff can be
        # modified
        all_subreactions = self._model.subreaction_data

        for codon, count in self.translation_data.codon_count.items():
            codon = codon.replace('U', 'T')
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

        # TODO: how many h/h2o molecules are exchanged when making peptide bond
        # tRNA charging requires 2 ATP per amino acid. TODO: tRNA demand
        # Translocation and elongation each use one GTP per event
        # (len(protein) - 1). Initiation and termination also each need
        # one GTP each. Therefore, the GTP hydrolysis resulting from
        # translation is 2 * len(protein)

        # tRNA + GTP -> tRNA_GTP
        # TODO: audit the numbers below. Check with subreaction data
        # TODO: handle in subreactions
        #new_stoichiometry["h2o_c"] -= 3 * protein_length
        #new_stoichiometry["h_c"] += 3 * protein_length
        #new_stoichiometry["pi_c"] += 3 * protein_length
        #new_stoichiometry["gtp_c"] -= 2 * protein_length
        #new_stoichiometry["gdp_c"] += 2 * protein_length
        #new_stoichiometry["atp_c"] -= 1 * protein_length
        #new_stoichiometry["adp_c"] += 1 * protein_length

        last_codon = self.translation_data.last_codon
        term_enzyme = model.translation_info["translation_terminators"].get(last_codon)
        if term_enzyme is not None:
            termination_subreaction_id = last_codon + '_' + term_enzyme + \
                '_mediated_termination'
            try:
                term_subreaction_data = \
                    all_subreactions.get_by_id(termination_subreaction_id)
            except KeyError as e:
                if verbose:
                    print("Termination subreaction '%s' for '%s' not found" %
                          (termination_subreaction_id, self.id))

            else:
                new_stoichiometry[term_subreaction_data.enzyme] -= \
                    mu / term_subreaction_data.keff / 3600.

        object_stoichiometry = self.get_components_from_ids(new_stoichiometry,
                                                            verbose=verbose)

        self.add_metabolites(object_stoichiometry,
                             combine=False, add_to_container_model=False)

        # add to biomass
        protein_mass = protein.formula_weight / 1000.  # kDa
        self.add_metabolites({self._model._protein_biomass: protein_mass},
                             combine=False, add_to_container_model=False)
        # RNA biomass
        RNA_mass = transcript.formula_weight / 1000.  # kDa
        self.add_metabolites(
            {self._model._RNA_biomass:
                 RNA_mass * (1 - deg_fraction) * RNA_amount},
            combine=False, add_to_container_model=False)


class tRNAChargingReaction(MEReaction):

    tRNAData = None

    def update(self, verbose=True):
        stoic = defaultdict(int)
        data = self.tRNAData

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

        stoic = self.add_modifications(self.tRNAData.id, stoic, tRNA_amount)

        object_stoichiometry = self.get_components_from_ids(stoic,
                                                            verbose=verbose)

        self.add_metabolites(object_stoichiometry, combine=False,
                             add_to_container_model=False)


class SummaryVariable(MEReaction):
    pass
