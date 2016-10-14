from __future__ import division

from warnings import warn
from collections import defaultdict, Counter
from itertools import product

from six import iteritems

from cobra import Reaction
from cobra.core.Formula import Formula

from cobrame.util.mass import *
from cobrame.util import mu, dogma
from sympy import Basic
from cobrame.core.Components import *


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
           Some processes (ie. tRNA charging) are reformulated such that other
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
                        mu / modification.keff / 3600. * scale * abs(count)
            elif type(modification.enzyme) == str:
                stoichiometry[modification.enzyme] -= \
                    mu / modification.keff / 3600. * scale * abs(count)

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
                    stoichiometry[enzyme] -= mu / subreaction_data.keff / \
                                             3600. * count
            elif type(subreaction_data.enzyme) == str:
                stoichiometry[subreaction_data.enzyme] -= \
                    mu / subreaction_data.keff / 3600. * count

            for met, stoich in subreaction_data.stoichiometry.items():
                stoichiometry[met] += count * stoich

        return stoichiometry

    def add_translocation_pathways(self, process_data_id, protein_id,
                                   stoichiometry):

        all_translocation = self._model.translocation_data
        process_info = self._model.process_data.get_by_id(process_data_id)
        protein = self._model.metabolites.get_by_id(protein_id)
        protein_length = len(protein.amino_acid_sequence)
        # Requirement of some translocation complexes vary depending
        # on protein being translocated
        multiplier_dict = self._model.global_info["translocation_multipliers"]

        for translocation, count in iteritems(process_info.translocation):
            translocation_data = all_translocation.get_by_id(translocation)
            for metabolite, amount in iteritems(translocation_data.stoichiometry):
                if translocation_data.length_dependent_energy:
                    stoichiometry[metabolite] += amount * count * \
                                                 protein_length
                else:
                    stoichiometry[metabolite] += amount * count

            for enzyme, enzyme_info in iteritems(
                    translocation_data.enzyme_dict):
                length_dependent = enzyme_info['length_dependent']
                fixed_keff = enzyme_info['fixed_keff']
                bnum = protein_id.replace('protein_', '')

                if multiplier_dict:
                    multiplier = multiplier_dict[enzyme].get(bnum, 1)
                else:
                    multiplier = 1

                if not length_dependent:
                    protein_length = 1.

                # TODO is this backward?
                # keff = translocation_data.keff
                if fixed_keff:
                    keff = 65.
                else:
                    keff = translocation_data.keff / protein_length

                enzyme_stoichiometry = multiplier * mu / keff / 3600. * count
                stoichiometry[enzyme] -= enzyme_stoichiometry

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
                    print("Created %s in %s" % (repr(new_met), repr(self)))
                object_stoichiometry[new_met] = value
                self._model.add_metabolites([new_met])
        return object_stoichiometry

    def check_ME_mass_balance(self):
        """Compute mass and charge balance for the reaction

        returns a dict of {element: amount} for unbalanced elements.
        "charge" is treated as an element in this dict
        This will be empty for balanced reactions.
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

                if met.id == 'generic_tRNA_GAG_glu__L_c':
                    reaction_element_dict['H'] -= coefficient
                    reaction_element_dict['O'] -= coefficient

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
    """Metabolic reaction including required enzymatic complex

    This reaction class's update function processes the a information contained
    in the complex data for the enzyme that catalyzes this reaction as well as
    the stoichiometric data which contains the stoichiometry of the metabolic
    conversion being performed (i.e. the stoichiometry of the M-model reaction
    analog)

    :param float keff:
        The keff couples enzymatic dilution to metabolic flux
    :param Boolean reverse:
        If True, the reaction corresponds to the reverse direction of the
        reaction. This is necessary since all reversible enzymatic reactions
        in an ME-model are broken into two irreversible reactions
    :param set complex_dilution_set:
        If an enzyme involved of this reaction acts as a "carrier" (i.e. enzyme
        that transfers a sidechain to another enzyme or metabolite)
    """

    @property
    def complex_data(self):
        return self._complex_data

    @complex_data.setter
    def complex_data(self, process_data):
        self._complex_data = process_data
        if process_data is not None:
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
    complex_dilution_set = set()

    def update(self, verbose=True):
        self.clear_metabolites()
        new_stoichiometry = defaultdict(float)
        stoichiometric_data = self.stoichiometric_data

        # Add complex if enzyme catalyzed
        if self.complex_data:
            new_stoichiometry[self.complex_data.complex.id] = \
                -mu / self.keff / 3600.  # s-1 / (3600 s/h)

        # Update new stoichiometry values
        sign = -1 if self.reverse else 1
        for component, value in iteritems(stoichiometric_data.stoichiometry):
            new_stoichiometry[component] += value * sign
            if component in self.complex_dilution_set:
                new_stoichiometry[component] += - mu / self.keff / 3600.

        # Convert component ids to cobra metabolites
        object_stoichiometry = self.get_components_from_ids(new_stoichiometry,
                                                            verbose=verbose)

        # Replace old stoichiometry with new one
        self.add_metabolites(object_stoichiometry,
                             add_to_container_model=False)

        # Set the bounds
        if self.reverse:
            self.lower_bound = max(0, -self.stoichiometric_data.upper_bound)
            self.upper_bound = max(0, -self.stoichiometric_data.lower_bound)
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
        self.clear_metabolites()
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
                             add_to_container_model=False)


class PostTranslationReaction(MEReaction):
    """
    Reaction class that includes all posttranslational modification reactions
    (translocation, protein folding, etc)

    There are often multiple different reactions/enzymes that can accomplish
    the same modification/function. In order to account for these and
    maintain one translation reaction per protein, these processes need to be
    modeled as separate reactions.

    Note: for now folding is handled as part of translation, but this will
    be updated once multiple folding pathways are accounted for.
    """

    _posttranslation_data = None

    @property
    def posttranslation_data(self):
        return self._posttranslation_data

    @posttranslation_data.setter
    def posttranslation_data(self, process_data):
        self._posttranslation_data = process_data
        process_data._parent_reactions.add(self.id)

    def update(self, verbose=True):
        self.clear_metabolites()
        stoichiometry = defaultdict(float)
        metabolites = self._model.metabolites
        posttranslation_data = self.posttranslation_data
        unprocessed_protein = posttranslation_data.unprocessed_protein_id
        processed_protein = posttranslation_data.processed_protein_id

        # folding properties
        folding_mechanism = posttranslation_data.folding_mechanism
        aggregation_propensity = posttranslation_data.aggregation_propensity
        scaling = posttranslation_data.propensity_scaling
        if folding_mechanism:
            T = str(self._model.global_info['temperature'])
            keq_folding = posttranslation_data.keq_folding[T]
            k_folding = posttranslation_data.k_folding[T] * 3600.  # in hr-1

        try:
            protein_met = metabolites.get_by_id(processed_protein)
        except KeyError:
            protein_met = ProcessedProtein(processed_protein,
                                           unprocessed_protein)
            self._model.add_metabolites(protein_met)

        stoichiometry = self.add_modifications(posttranslation_data.id,
                                               stoichiometry)
        stoichiometry = self.add_subreactions(posttranslation_data.id,
                                              stoichiometry)

        if posttranslation_data.translocation:
            stoichiometry = self.add_translocation_pathways(posttranslation_data.id,
                                                            unprocessed_protein,
                                                            stoichiometry)
        if folding_mechanism == 'folding_spontaneous':
            dilution = (keq_folding + mu / k_folding)
            stoichiometry[unprocessed_protein] -= (dilution + 1.)
            stoichiometry[protein_met.id] += 1.

        elif folding_mechanism:
            dilution = aggregation_propensity * scaling * (keq_folding + 1.)+1.
            stoichiometry[unprocessed_protein] -= (1. / dilution + 1.)

            stoichiometry[protein_met.id] += 1. / dilution
            stoichiometry[protein_met.id.replace('_folded', '')] += (1.)
        else:
            stoichiometry[unprocessed_protein] = -1.
            stoichiometry[protein_met.id] = 1.

        # Add surface area constraints for all translocated protein
        surface_area = posttranslation_data.surface_area
        if surface_area:
            for SA, value in surface_area.items():
                try:
                    SA_constraint = metabolites.get_by_id(SA)
                except KeyError:
                    warn('Constraint %s added to model' % SA)
                    SA_constraint = Constraint(SA)
                    self._model.add_metabolites([SA_constraint])

                stoichiometry[SA_constraint.id] += value

        object_stoichiometry = self.get_components_from_ids(stoichiometry,
                                                            verbose=verbose)

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
        self.clear_metabolites()
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
            k_RNAP = (mu * c_ribo * kt / (mu + kt * r0)) * 3  # (3*k_ribo) hr-1
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

        # Add modifications and subreactions to reaction stoichiometry
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
        tRNA_mass = 0.
        rRNA_mass = 0.
        ncRNA_mass = 0.
        mRNA_mass = 0.
        for met, v in new_stoich.items():
            if v < 0 or not hasattr(met, "RNA_type"):
                continue
            if met.RNA_type == 'tRNA':
                tRNA_mass += met.mass  # kDa
            if met.RNA_type == 'rRNA':
                rRNA_mass += met.mass  # kDa
            if met.RNA_type == 'ncRNA':
                ncRNA_mass += met.mass  # kDa
            if met.RNA_type == 'mRNA':
                mRNA_mass += met.mass  # kDa

        if tRNA_mass > 0:
            self.add_metabolites({self._model._tRNA_biomass: tRNA_mass},
                                 combine=False, add_to_container_model=False)
        if rRNA_mass > 0:
            self.add_metabolites({self._model._rRNA_biomass: rRNA_mass},
                                 combine=False, add_to_container_model=False)
        if ncRNA_mass > 0:
            self.add_metabolites({self._model._ncRNA_biomass: ncRNA_mass},
                                 combine=False, add_to_container_model=False)
        if mRNA_mass > 0:
            self.add_metabolites({self._model._mRNA_biomass: mRNA_mass},
                                 combine=False, add_to_container_model=False)


class GenericFormationReaction(MEReaction):
    pass


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
        self.clear_metabolites()
        protein_id = self.translation_data.protein
        mRNA_id = self.translation_data.mRNA
        protein_length = len(self.translation_data.amino_acid_sequence)
        model = self._model
        metabolites = self._model.metabolites
        new_stoichiometry = defaultdict(int)
        elements = defaultdict(int)

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

        # -----------------Add Amino Acids----------------------------------
        for aa, value in self._translation_data.amino_acid_count.items():
            new_stoichiometry[aa] = -value

        # -----------------Add Ribosome Coupling----------------------------
        try:
            ribosome = metabolites.get_by_id("ribosome")
        except KeyError:
            warn("ribosome not found")
        else:
            k_ribo = mu * c_ribo * kt / (mu + kt * r0)  # in hr-1
            coupling = -protein_length * mu / k_ribo
            new_stoichiometry[ribosome.id] = coupling

        # -------------------Add mRNA Coupling------------------------------
        try:
            transcript = metabolites.get_by_id(mRNA_id)
        except KeyError:
            # If transcript not found add to the model as the mRNA_id
            warn("transcript '%s' not found" % mRNA_id)
            transcript = TranscribedGene(mRNA_id)
            model.add_metabolites(transcript)

        # Calculate coupling constraints for mRNA and degradation
        k_mRNA = mu * c_mRNA * kt / (mu + kt * r0) * 3.  # in hr-1 TODO should be *3
        RNA_amount = mu / k_mRNA
        #deg_fraction = 3. * k_deg / (3. * k_deg + mu)
        deg_fraction = 1. / (k_deg) if k_deg != 0 else 0
        deg_amount = deg_fraction * RNA_amount

        # Add mRNA coupling to stoichiometry
        new_stoichiometry[transcript.id] = -RNA_amount

        # ---------------Add Degradation Requirements -------------------------
        # Add degraded nucleotides to stoichiometry
        for nucleotide, count in transcript.nucleotide_count.items():
            new_stoichiometry[nucleotide] += count * deg_amount

        # ATP hydrolysis required for cleaving
        nucleotide_length = len(transcript.nucleotide_sequence)
        hydrolysis_amount = (nucleotide_length - 1) / 4. * deg_amount
        atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1,
                          'pi_c': 1, 'h_c': 1}
        for metabolite, value in atp_hydrolysis.items():
            new_stoichiometry[metabolite] = hydrolysis_amount * value

        # Add degradosome coupling, if known
        try:
            RNA_degradosome = metabolites.get_by_id("RNA_degradosome")
        except KeyError:
            warn("RNA_degradosome not found")
        else:
            deg_coupling = -deg_amount * mu / 65. / 3600. # keff of degradosome
            new_stoichiometry[RNA_degradosome.id] = deg_coupling

        # --------------- Add Protein to Stoichiometry ------------------------
        # Added protein to model if not already included
        try:
            protein = metabolites.get_by_id(protein_id)
        except KeyError:
            protein = TranslatedGene(protein_id)
            model.add_metabolites(protein)
        new_stoichiometry[protein.id] = 1

        # ------------------ Elongation Reactions------------------------
        # We add in a "generic tRNA" for each amino acid. The production
        # of this generic tRNA represents the production of enough of any tRNA
        # (according to its keff) for that amino acid to catalyze addition
        # of a single amino acid to a growing peptide.

        # Add the subreaction_data associated with each tRNA/AA addition
        # The subreaction data itself is added in building.py

        # Addition of translation start subreaction handled in buidling.py

        # Correction for formylmethionine start codon
        elements["C"] += 1
        elements["O"] += 1

        # tRNA charging requires 2 ATP per amino acid.
        # Translocation and elongation each use one GTP per event
        # (len(protein) - 1). Initiation and termination also each need
        # one GTP each. Therefore, the GTP hydrolysis resulting from
        # translation is 2 * len(protein)
        # tRNA + GTP -> tRNA_GTP

        # Add transcription termination as subreaction, handled in building.py

        # ------- Convert ids to metabolites and add to model -----------------
        # add subreactions to stoichiometry
        new_stoichiometry = self.add_subreactions(self.translation_data.id,
                                                  new_stoichiometry)
        new_stoichiometry = self.add_modifications(self.translation_data.id,
                                                   new_stoichiometry)
        # convert metabolite ids to cobra metabolites
        object_stoichiometry = self.get_components_from_ids(new_stoichiometry,
                                                            verbose=verbose)
        # add metabolites to reaction
        self.add_metabolites(object_stoichiometry,
                             combine=False, add_to_container_model=False)

        # ------------------ Add biomass constraints --------------------------
        # add to biomass
        protein_mass = protein.mass  # kDa
        self.add_metabolites({self._model._protein_biomass: protein_mass},
                             combine=False, add_to_container_model=False)
        # RNA biomass
        mRNA_mass = transcript.mass  # kDa
        self.add_metabolites(
            {self._model._mRNA_biomass: (-mRNA_mass * deg_amount)},
            combine=False, add_to_container_model=False)

        # -------------Update Element Dictionary and Formula-------------------
        aa_count = self.translation_data.amino_acid_count
        for aa_name, value in iteritems(aa_count):
            for e, n in iteritems(metabolites.get_by_id(aa_name).elements):
                elements[e] += n * value
        # subtract water from composition
        elements["H"] -= (protein_length - 1) * 2
        elements["O"] -= protein_length - 1
        protein.formula = "".join(stringify(e, n)
                                  for e, n in sorted(iteritems(elements)))


class tRNAChargingReaction(MEReaction):

    tRNAData = None

    def update(self, verbose=True):
        self.clear_metabolites()
        new_stoichiometry = defaultdict(float)
        data = self.tRNAData

        # set tRNA coupling parameters
        m_tRNA = self._model.global_info['m_tRNA']
        m_aa = self._model.global_info['m_aa']
        f_tRNA = self._model.global_info['f_tRNA']
        kt = self._model.global_info['kt']  # hr-1
        r0 = self._model.global_info['r0']
        c_tRNA = m_tRNA / m_aa / f_tRNA

        # The meaning of a generic tRNA is described in the
        # TranslationReaction comments
        generic_tRNA = "generic_tRNA_" + data.codon + "_" + data.amino_acid
        new_stoichiometry[generic_tRNA] = 1

        # Compute tRNA (and amino acid) coupling and add to stoichiometry
        tRNA_keff = c_tRNA * kt * mu / (mu + r0 * kt)  # per hr
        tRNA_amount = mu / tRNA_keff
        new_stoichiometry[data.RNA] = -tRNA_amount
        new_stoichiometry[data.amino_acid] = -tRNA_amount

        # Add synthetase coupling and enzyme, if known
        synthetase_amount = mu / data.synthetase_keff / \
                            3600. * (1 + tRNA_amount)
        if data.synthetase is not None:
            new_stoichiometry[data.synthetase] = -synthetase_amount

        # Add tRNA modifications to stoichiometry
        new_stoichiometry = self.add_modifications(self.tRNAData.id,
                                                   new_stoichiometry,
                                                   tRNA_amount)

        # Convert component ids to cobra metabolites
        object_stoichiometry = self.get_components_from_ids(new_stoichiometry,
                                                            verbose=verbose)

        # Replace reaction stoichiometry with updated stoichiometry
        self.add_metabolites(object_stoichiometry,
                             add_to_container_model=False)


class SummaryVariable(MEReaction):
    pass
