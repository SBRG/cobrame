from __future__ import division

from warnings import warn
from collections import defaultdict, Counter
from itertools import product

from six import iteritems

from cobra import Reaction
from cobra.core.Formula import Formula

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
                        mu / modification.keff / 3600. * scale
            elif type(modification.enzyme) == str:
                stoichiometry[modification.enzyme] -= \
                    mu / modification.keff / 3600. * scale

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
        metabolites = self._model.metabolites
        coefficient_dict = {}
        reaction_element_dict = defaultdict(int)
        for metabolite, coefficient in self._metabolites.items():
            if metabolite.id == 'biomass':
                continue
            if isinstance(coefficient, Basic):
                if isinstance(metabolite, Ribosome):
                    coefficient = 0
                elif isinstance(metabolite, RNAP):
                    coefficient = 0
                elif isinstance(metabolite, TranscribedGene) and isinstance(self, tRNAChargingReaction):
                    coefficient = 0
                else:
                    coefficient = coefficient.subs(mu, 0)
            coefficient_dict[metabolite] = coefficient
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += \
                    coefficient * metabolite.charge
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coefficient * amount

        # Don't account for remaining nucleotides on TU if they are excised
        # Only include intergenic bases for solely mRNA coding TUs.
        # TUs including both mRNA and t(r)RNA are already mass balanced.
        if isinstance(self, TranscriptionReaction):
            data_id = self.id.replace('transcription_', '')
            TU = self._model.transcription_data.get_by_id(data_id)
            if len(TU.excised_bases) == 0:
                TU_seq = TU.nucleotide_sequence
                TU_seq_red = TU_seq
                left_pos_list = []
                right_pos_list = []
                for RNA in TU.RNA_products:
                    met = metabolites.get_by_id(RNA)
                    strand = met.strand
                    dna_seq = met.nucleotide_sequence
                    function_3 = lambda x, y: dna_seq[x+1:-y] \
                        if strand == '+' else dna_seq[y:-(x+1)]
                    function_1 = lambda x: dna_seq[x:] \
                        if strand == '+' else dna_seq[:-(x+1)]
                    function_2 = lambda x: dna_seq[:-x] \
                        if strand == '+' else dna_seq[x+1:]
                    if dna_seq in TU_seq and dna_seq not in TU_seq_red:
                        for left, right in zip(left_pos_list, right_pos_list):
                            if met.left_pos < right < met.right_pos:
                                change = right - met.left_pos
                                new_seq = function_1(change)
                                TU_seq_red = TU_seq_red.replace(new_seq, '')
                            elif met.left_pos < left < met.right_pos:
                                change = met.right_pos - left
                                new_seq = function_2(change)
                                TU_seq_red = TU_seq_red.replace(new_seq, '')
                        for left, right in product(left_pos_list, right_pos_list):
                            if met.left_pos < right and met.right_pos > left:
                                left_change = right - met.left_pos
                                right_change = met.right_pos - left
                                new_seq = function_3(left_change, right_change)
                                TU_seq_red = TU_seq_red.replace(new_seq, '')
                    else:
                        TU_seq_red = TU_seq_red.replace(dna_seq, '')
                    left_pos_list.append(met.left_pos)
                    right_pos_list.append(met.right_pos)
                for nucleotide, value in iteritems(Counter(TU_seq_red)):
                    nucleotide = nucleotide.replace('T', 'U')
                    met = metabolites.get_by_id(nucleotide.lower() + 'mp_c')
                    for element, amount in iteritems(met.elements):
                        reaction_element_dict[element] += value * amount

        if isinstance(self, MetabolicReaction) or isinstance(self, tRNAChargingReaction):
            mod_formula_dict = self._model.global_info['modification_formulas']
            for met in self.metabolites:
                coefficient = coefficient_dict[met]
                if not isinstance(met, Complex) or coefficient == 0:
                    continue

                mods = met.id.split('_mod_')
                if len(mods) == 1:
                    continue

                start_mod = 1

                # elements from complexes with formulas are handled above
                if not met.formula:
                    met_list = met.id.split('_mod_')
                    cplx_1 = met_list[0] + '_mod_' + met_list[1]

                    try:
                        complex_met = metabolites.get_by_id(cplx_1)
                        formula = complex_met.formula
                        start_mod = 2
                    except:
                        formula = None

                    if not formula:
                        start_mod = 1
                        complex_met = metabolites.get_by_id(met_list[0])

                    for element, amount in iteritems(complex_met.elements):
                        reaction_element_dict[element] += coefficient * amount

                    for mod in mods[start_mod:]:
                        if mod == 'Oxidized':
                            pass
                            #reaction_element_dict['H'] += 2. * coefficient
                        if len(mod.split(':')) == 1:
                            value = 1
                            mod_met = mod
                        else:
                            value, mod_met = mod.split(':')
                        try:
                            formula = mod_formula_dict[mod_met]['formula']
                            formula_obj = Formula(str(formula))
                            elements = formula_obj.elements
                        except KeyError:
                            elements = metabolites.get_by_id(mod_met + '_c').elements
                        for element, amount in iteritems(elements):
                            reaction_element_dict[element] += coefficient * amount * float(value)


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
        self._metabolites.clear()
        metabolites = self._model.metabolites
        new_stoichiometry = defaultdict(float)
        if self.complex_data:
            new_stoichiometry[self.complex_data.complex] = \
                -mu / self.keff / 3600.  # s-1 / (3600 s/h)
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
        self._metabolites.clear()
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
        self._metabolites.clear()
        TU_id = self.transcription_data.id
        stoichiometry = defaultdict(int)
        TU_length = len(self.transcription_data.nucleotide_sequence)
        RNA_polymerase = self.transcription_data.RNA_polymerase
        metabolites = self._model.metabolites
        # Set Parameters
        kt = self._model.global_info['kt']
        r0 = self._model.global_info['r0']
        m_rr = self._model.global_info['m_rr']
        f_rRNA = self._model.global_info['f_rRNA']
        m_aa = self._model.global_info['m_aa']
        c_ribo = m_rr / f_rRNA / m_aa

        try:
            RNAP = self._model.metabolites.get_by_id(RNA_polymerase)
        except KeyError:
            warn("RNA Polymerase (%s) not found" % RNA_polymerase)
        else:
            k_RNAP = (mu * c_ribo * kt / (mu + kt * r0)) * 3  # in hr-1 3*k_ribo
            coupling = -TU_length * mu / k_RNAP
            stoichiometry[RNAP.id] = coupling

        # This is necessary as opposed to using get_components_from_ids in
        # so that self.excised_bases can be used for "dummy_RNA"
        # TODO Can we fix this?
        for transcript_id in self.transcription_data.RNA_products:
            if transcript_id not in metabolites:
                warn("Transcript (%s) added to model " % transcript_id)
                transcript = TranscribedGene(transcript_id)
                self._model.add_metabolites([transcript])
                demand_reaction = Reaction("DM_" + transcript_id)
                self._model.add_reaction(demand_reaction)
                demand_reaction.add_metabolites({transcript_id: -1})
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

        stoichiometry["ppi_c"] += TU_length
        stoichiometry["h2o_c"] -= TU_length
        stoichiometry["h_c"] += TU_length
        # 5' had a triphosphate, but this is removed when excising
        # Is this universally true?
        if sum(self.transcription_data.excised_bases.values()) > 0:
            stoichiometry["ppi_c"] += 1

        new_stoich = \
            self.get_components_from_ids(stoichiometry, verbose=verbose,
                                         default_type=TranscribedGene)

        self.add_metabolites(new_stoich, combine=False,
                             add_to_container_model=False)

        # mRNA biomass contribution handled in translation
        tRNA_mass = 0
        rRNA_mass = 0
        ncRNA_mass = 0
        for met, v in new_stoich.items():
            if v < 0 or not hasattr(met, "RNA_type"):
                continue
            if met.RNA_type == 'tRNA':
                tRNA_mass += met.mass  # kDa
            if met.RNA_type == 'rRNA':
                rRNA_mass += met.mass  # kDa
            if met.RNA_type == 'ncRNA':
                ncRNA_mass += met.mass  # kDa
        if tRNA_mass > 0:
            self.add_metabolites({self._model._tRNA_biomass: tRNA_mass},
                                 combine=False, add_to_container_model=False)
        if rRNA_mass > 0:
            self.add_metabolites({self._model._rRNA_biomass: rRNA_mass},
                                 combine=False, add_to_container_model=False)
        if ncRNA_mass > 0:
            self.add_metabolites({self._model._ncRNA_biomass: ncRNA_mass},
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
        self._metabolites.clear()
        protein_id = self.translation_data.protein
        mRNA_id = self.translation_data.mRNA
        protein_length = len(self.translation_data.amino_acid_sequence)
        model = self._model
        metabolites = self._model.metabolites
        new_stoichiometry = defaultdict(int)
        # Set Parameters
        kt = self._model.global_info['kt']
        k_deg = self._model.global_info['k_deg']
        r0 = self._model.global_info['r0']
        m_rr = self._model.global_info['m_rr']
        f_rRNA = self._model.global_info['f_rRNA']
        m_aa = self._model.global_info['m_aa']

        m_nt = self._model.global_info['m_nt']
        f_mRNA = self._model.global_info['f_mRNA']

        c_ribo = m_rr / f_rRNA / m_aa
        c_mRNA = m_nt / f_mRNA / m_aa

        try:
            ribosome = metabolites.get_by_id("ribosome")
        except KeyError:
            warn("ribosome not found")
        else:
            k_ribo = mu * c_ribo * kt / (mu + kt * r0)  # in hr-1
            coupling = -protein_length * mu / k_ribo
            new_stoichiometry[ribosome.id] = coupling

        # Not all genes have annotated TUs. For these add TU as the bnumber
        try:
            transcript = metabolites.get_by_id(mRNA_id)
        except KeyError:
            warn("transcript '%s' not found" % mRNA_id)
            transcript = TranscribedGene(mRNA_id)
            model.add_metabolites(transcript)

        k_mRNA = mu * c_mRNA * kt / (mu + kt * r0)  # in hr-1
        RNA_amount = mu / k_mRNA

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

        # TODO account of selenocysteine in formula

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
                if abbreviated_aa == "*":
                    break
                aa = dogma.amino_acids[abbreviated_aa].split('_')[0]
            codon = codon.replace('T', 'U')
            subreaction_id = aa + '_addition_at_' + codon
            try:
                self.translation_data.subreactions[subreaction_id] = count
            except KeyError:
                warn('subreaction %s not in model' % subreaction_id)

        protein.formula = "".join(stringify(e, n)
                                  for e, n in sorted(iteritems(elements)))

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
        term_enzyme = model.global_info["translation_terminators"].get(last_codon)
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
        protein_mass = protein.mass  # kDa
        self.add_metabolites({self._model._protein_biomass: protein_mass},
                             combine=False, add_to_container_model=False)
        # RNA biomass
        mRNA_mass = transcript.mass  # kDa
        self.add_metabolites(
            {self._model._mRNA_biomass:
                 mRNA_mass * (1 - deg_fraction) * RNA_amount},
            combine=False, add_to_container_model=False)


class tRNAChargingReaction(MEReaction):

    tRNAData = None

    def update(self, verbose=True):
        self._metabolites.clear()
        stoic = defaultdict(int)
        data = self.tRNAData

        # set tRNA coupling parameters
        m_tRNA = self._model.global_info['m_tRNA']
        m_aa = self._model.global_info['m_aa']
        f_tRNA = self._model.global_info['f_tRNA']
        kt = self._model.global_info['kt']  # hr-1
        r0 = self._model.global_info['r0']
        c_tRNA = m_tRNA / m_aa / f_tRNA
        tRNA_keff = c_tRNA * kt * mu / (mu + r0 * kt)  # per hr

        # If the generic tRNA does not exist, create it now. The meaning of
        # a generic tRNA is described in the TranslationReaction comments
        generic_tRNA = "generic_tRNA_" + data.codon + "_" + data.amino_acid

        stoic[generic_tRNA] = 1

        # compute what needs to go into production of a generic tRNA
        tRNA_amount = mu / tRNA_keff
        synthetase_amount = mu / data.synthetase_keff / \
                            3600. * (1 + tRNA_amount)
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
