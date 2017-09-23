from collections import defaultdict

from six import iteritems
import cobra
import numpy as np

from cobrame.core.MEReactions import *
from cobrame.core.Components import GenerictRNA
from cobrame.util.dogma import *
from cobrame.util.massbalance import elements_to_formula


class ProcessData(object):
    """Generic class for holding information about a process

    ME reactions are built from information in these objects

    """

    def __init__(self, id, model):
        self.id = id
        self._model = model
        # parents need to be updated every time the process is updated
        # a parent must have an update method
        self._parent_reactions = set()
        model.process_data.append(self)

    @property
    def model(self):
        return self._model

    @property
    def parent_reactions(self):
        return {self._model.reactions.get_by_id(i)
                for i in self._parent_reactions}

    def _update_parent_reactions(self):
        reactions = self._model.reactions
        for i in self._parent_reactions:
            reactions.get_by_id(i).update()

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))


class StoichiometricData(ProcessData):
    """Encodes the stoichiometry for a metabolic reaction.

    This stoichiometry and ID will often link directly to M-model reactions

    :param dict _stoichiometry:
        Dictionary of {metabolite_id: stoichiometry}

    :param defaultdict subreactions:
        This should rarely contain anything. Only in special cases where
        a second diluted enzyme needs added to metabolic reaction
    """
    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.stoichiometric_data.append(self)
        self._stoichiometry = {}
        self.subreactions = defaultdict(int)
        self.lower_bound = 0.
        self.upper_bound = 1000.

    @property
    def stoichiometry(self):
        return self._stoichiometry


class ModificationData(ProcessData):
    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.modification_data.append(self)
        self.stoichiometry = {}
        self.enzyme = None
        self.keff = 65.
        self._element_contribution = {}

    @property
    def element_contribution(self):
        if self._element_contribution:
            return self._element_contribution

        # Return "trivial" cases (only one modifying metabolite in the
        # reactants and no products) without warning
        elif len(self.stoichiometry) == 1 and \
                        list(self.stoichiometry.values())[0] < 0:
            return self.calculate_element_contribution()
        elif self.calculate_element_contribution():
            warn('No element contribution input for modification (%s), '
                 'calculating based on stoichiometry instead' % self.id)
            return self.calculate_element_contribution()
        else:
            return {}

    def calculate_element_contribution(self):
        elements = defaultdict(int)
        for met, coefficient in iteritems(self.stoichiometry):
            met_obj = self._model.metabolites.get_by_id(met)
            # elements lost in conversion are added to complex, protein, etc.
            if not met_obj.elements and not isinstance(met_obj, GenerictRNA):
                warn('met (%s) does not have formula' % met_obj.id)
            for e, n in iteritems(met_obj.elements):
                elements[e] -= n * coefficient
        return elements

    def calculate_biomass_contribution(self):
        elements = self.calculate_element_contribution()

        # Create temporary metabolite for calculating formula weight
        tmp_met = cobra.Metabolite('mass')
        elements_to_formula(tmp_met, elements)

        return tmp_met.formula_weight

    def get_complex_data(self):
        for i in self._model.complex_data:
            if self.id in i.modifications:
                yield i


class SubreactionData(ProcessData):
    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.subreaction_data.append(self)
        self.stoichiometry = {}
        self.enzyme = None
        self.keff = 65.
        self._element_contribution = {}

    @property
    def element_contribution(self):
        if self._element_contribution:
            return self._element_contribution
        else:
            contribution = {i: v for i, v in
                            iteritems(self.calculate_element_contribution())
                            if v}

        # Return "trivial" cases (only one modifying metabolite in the
        # reactants and no products) without warning
        if len(self.stoichiometry) == 1 and \
                        list(self.stoichiometry.values())[0] < 0:
            return self.calculate_element_contribution()
        elif contribution:
            warn('No element contribution input for subreaction (%s), '
                 'calculating based on stoichiometry instead' % self.id)
            return self.calculate_element_contribution()
        else:
            return {}

    def calculate_element_contribution(self):
        elements = defaultdict(int)
        for met, coefficient in iteritems(self.stoichiometry):
            met_obj = self._model.metabolites.get_by_id(met)
            # elements lost in conversion are added to complex, protein, etc.
            if not met_obj.elements and not isinstance(met_obj, GenerictRNA):
                warn('met (%s) does not have formula' % met_obj.id)
            for e, n in iteritems(met_obj.elements):
                elements[e] -= n * coefficient
        return elements

    def calculate_biomass_contribution(self):
        elements = self.calculate_element_contribution()

        # Create temporary metabolite for calculating formula weight
        tmp_met = cobra.Metabolite('mass')
        elements_to_formula(tmp_met, elements)

        return tmp_met.formula_weight

    def get_complex_data(self):
        for i in self._model.complex_data:
            if self.id in i.subreactions:
                yield i


class ComplexData(ProcessData):

    """Contains all information associated with the formation of an
    enzyme complex

    :param dict stoichiometry:
        Dictionary containing {protein_id: count} for all protein subunits
        comprising enzyme complex

    :param dict modifications:
        Dictionary of {modfication_data_id: count} for all protein
        modifications. This can include cofactor/prosthetic group binding
        or enzyme side group addition.

    """

    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.complex_data.append(self)
        # {Component.id: stoichiometry}
        self.stoichiometry = defaultdict(float)
        self.chaperones = {}
        # {ModificationData.id : number}
        self.modifications = {}
        # Forming some metacomplexes occur in multiple steps
        self.subreactions = {}
        self._complex_id = None  # assumed to be the same as id if None

    @property
    def formation(self):
        """a read-only link to the formation reaction object"""
        try:
            return self._model.reactions.get_by_id("formation_" + self.id)
        except KeyError:
            return None

    @property
    def complex(self):
        """a read-only link to the complex object"""
        return self._model.metabolites.get_by_id(self.complex_id)

    @property
    def complex_id(self):
        """There are cases where multiple equivalent processes can result in
        the same final complex. This allows the equivalent final complex
        complex_id to be queried. This only needs set in the above case"""
        return self.id if self._complex_id is None else self._complex_id

    @complex_id.setter
    def complex_id(self, value):
        self._complex_id = None if value == self.id else value

    def create_complex_formation(self, verbose=True):
        """creates a complex formation reaction

        This assumes none exists already. Will create a reaction (prefixed by
        "formation") which forms the complex"""
        formation_id = "formation_" + self.id
        if formation_id in self._model.reactions:
            raise ValueError("reaction %s already in model" % formation_id)
        formation = ComplexFormation(formation_id)
        formation.complex_data_id = self.id
        formation._complex_id = self.complex_id
        self._model.add_reaction(formation)
        formation.update(verbose=verbose)


class TranscriptionData(ProcessData):
    def __init__(self, id, model, RNA_products=set()):
        ProcessData.__init__(self, id, model)
        model.transcription_data.append(self)
        self.nucleotide_sequence = ''
        self.RNA_products = RNA_products
        self.RNA_polymerase = ''
        # {ModificationData.id : number}
        self.modifications = defaultdict(int)
        self.subreactions = defaultdict(int)

    @property
    def nucleotide_count(self):
        return {transcription_table[i]: self.nucleotide_sequence.count(i)
                for i in ["A", "T", "G", "C"]}

    @property
    def RNA_types(self):
        return (self._model.metabolites.get_by_id(i).RNA_type for i in
                self.RNA_products)

    @property
    def mass(self):
        return compute_RNA_mass(self.nucleotide_sequence, self.excised_bases)

    # Bases to be excised due to splicing tRNA, rRNA, or ncRNA in RNA_products
    # i.e. {"amp_c": 10, "gmp_c": 11, "ump_c": 9, "cmp_c": 11}
    @property
    def excised_bases(self):
        # Skip if TU does not have any annotated RNA Products
        if len(self.RNA_products) == 0:
            return {}

        # Skip if TU only codes for mRNA
        RNA_types = set(self.RNA_types)
        if RNA_types == {"mRNA"}:
            return {}

        # Get dictionary of all nucleotide counts for TU
        seq = self.nucleotide_sequence
        counts = {i: seq.count(i) for i in ("A", "T", "G", "C")}

        # Subtract bases contained in RNA_product from dictionary
        metabolites = self._model.metabolites
        for product_id in self.RNA_products:
            gene_seq = metabolites.get_by_id(product_id).nucleotide_sequence
            for b in ("A", "T", "G", "C"):
                counts[b] -= gene_seq.count(b)

        # First base being a triphosphate will be handled by the reaction
        # producing an extra ppi during transcription. But generally, we add
        # triphosphate bases when transcribing, but excise monophosphate bases.
        monophosphate_counts = {dogma.transcription_table[k].replace("tp_c",
                                                                     "mp_c"): v
                                for k, v in iteritems(counts)}

        return monophosphate_counts

    @property
    def codes_stable_rna(self):
        # Return true if a stable RNA is in RNA_products
        has_stable_rna = False
        for RNA in self.RNA_products:
            try:
                gene = self._model.metabolites.get_by_id(RNA)
            except KeyError:
                pass
            else:
                if gene.RNA_type in ['tRNA', 'rRNA', 'ncRNA']:
                    has_stable_rna = True
        return has_stable_rna


class GenericData(ProcessData):
    def __init__(self, id, model, component_list):
        if not id.startswith("generic_"):
            warn("best practice for generic id to start with generic_")
        ProcessData.__init__(self, id, model)
        model.generic_data.append(self)
        self.component_list = component_list

    def create_reactions(self):
        model = self._model
        try:
            generic_metabolite = model.metabolites.get_by_id(self.id)
        except KeyError:
            generic_metabolite = GenericComponent(self.id)
            model.add_metabolites([generic_metabolite])
        for c_id in self.component_list:
            reaction_id = c_id + "_to_" + self.id
            try:
                reaction = model.reactions.get_by_id(reaction_id)
            except KeyError:
                reaction = GenericFormationReaction(reaction_id)
                model.add_reaction(reaction)
            stoic = {generic_metabolite: 1,
                     model.metabolites.get_by_id(c_id): -1}
            reaction.add_metabolites(stoic, combine=False)


class TranslationData(ProcessData):
    def __init__(self, id, model, mRNA, protein):
        ProcessData.__init__(self, id, model)
        model.translation_data.append(self)
        self.mRNA = mRNA
        self.protein = protein
        self.subreactions = defaultdict(int)
        self.modifications = defaultdict(int)
        self._amino_acid_sequence = ""
        self.nucleotide_sequence = ""
        self.term_enzyme = None

    @property
    def amino_acid_sequence(self):
        codons = (self.nucleotide_sequence[i: i + 3]
                  for i in range(0, (len(self.nucleotide_sequence)), 3))
        amino_acid_sequence = ''.join(codon_table.get(i, "K") for i in codons)
        amino_acid_sequence = amino_acid_sequence.rstrip("*")
        if not amino_acid_sequence.startswith('M'):
            start_codons = self._model.global_info.get("met_start_codons", {})
            if self.first_codon not in start_codons:
                warn("%s starts with '%s' which is not a start codon" %
                     (self.mRNA, self.first_codon))
            amino_acid_sequence = 'M' + ''.join(amino_acid_sequence[1:])
        if "*" in amino_acid_sequence:
            amino_acid_sequence = amino_acid_sequence.replace('*', 'C')
        return amino_acid_sequence

    @property
    def last_codon(self):
        return self.nucleotide_sequence[-3:].replace('T', 'U')

    @property
    def first_codon(self):
        return self.nucleotide_sequence[:3].replace('T', 'U')

    def _itercodons(self):
        None
        yield "AAU"

    @property
    def codon_count(self):
        # exclude the last three stop codons from count
        codons = (self.nucleotide_sequence[i: i + 3]
                  for i in range(0, (len(self.nucleotide_sequence)-3), 3))
        codon_count = defaultdict(int)
        for i in codons:
            codon_count[i.replace('T', 'U')] += 1

        # Remove one methionine (AUG) from codon count to account for start
        first_codon = self.first_codon
        start_codons = self._model.global_info.get("met_start_codons", {})
        if first_codon in start_codons:
            codon_count[first_codon] -= 1
        else:
            warn("%s starts with '%s' which is not a start codon" %
                 (self.mRNA, first_codon))
        return codon_count

    @property
    def amino_acid_count(self):
        """count of each amino acid in the protein"""
        aa_count = defaultdict(int)
        for i in self.amino_acid_sequence:
            aa_count[amino_acids[i]] += 1
        return aa_count

    @property
    def mass(self):
        """mass in kDa"""
        return compute_protein_mass(self.amino_acid_count)

    @property
    def elongation_subreactions(self):
        global_info = self._model.global_info
        elongation_subreactions = {}
        # Add the subreaction_data associated with each tRNA/AA addition
        for codon, count in iteritems(self.codon_count):
            codon = codon.replace('U', 'T')
            if codon == 'TGA':
                print('Adding selenocystein for %s' % self.id)
                aa = 'sec'
            else:
                abbreviated_aa = dogma.codon_table[codon]
                if abbreviated_aa == "*":
                    break
                # Filter out the compartment and stereochemistry from aa id
                aa = dogma.amino_acids[abbreviated_aa].split('_')[0]
            codon = codon.replace('T', 'U')
            subreaction_id = aa + '_addition_at_' + codon
            try:
                self._model.subreaction_data.get_by_id(subreaction_id)
            except KeyError:
                warn('elongation subreaction %s not in model' % subreaction_id)
            else:
                elongation_subreactions[subreaction_id] = count

        # Some additional enzymatic processes are required for each amino acid
        # addition during translation elongation
        translation_elongation_subreactions = \
            global_info.get('translation_elongation_subreactions', [])
        for subreaction_id in translation_elongation_subreactions:
            try:
                self._model.subreaction_data.get_by_id(subreaction_id)
            except KeyError:
                warn('elongation subreaction %s not in model' %
                     subreaction_id)
            else:
                # No elongation subreactions needed for start codon
                elongation_subreactions[subreaction_id] = \
                    len(self.amino_acid_sequence) - 1.

        return elongation_subreactions

    @property
    def translation_start_subreactions(self):
        # Read-only link to list of start subreactions for organism, if present
        # in model
        start_subreactions = []
        subreactions = \
            self._model.global_info.get('translation_start_subreactions', [])
        for subreaction_id in subreactions:
            try:
                self._model.subreaction_data.get_by_id(subreaction_id)
            except KeyError:
                warn('initiation subreaction %s not in model' %
                     subreaction_id)
            else:
                start_subreactions.append(subreaction_id)

        return start_subreactions

    @property
    def translation_termination_subreactions(self):
        termination_subreactions = []
        all_subreactions = self._model.subreaction_data
        global_info = self._model.global_info
        last_codon = self.last_codon
        term_enzyme = global_info["translation_terminators"].get(last_codon)
        if term_enzyme:
            termination_subreaction_id = last_codon + '_' + term_enzyme + \
                                         '_mediated_termination'
            try:
                all_subreactions.get_by_id(termination_subreaction_id)
            except KeyError:
                warn("Termination subreaction '%s' not in model" %
                     (termination_subreaction_id))
            else:
                termination_subreactions.append(termination_subreaction_id)
        return termination_subreactions


class tRNAData(ProcessData):
    synthetase = None
    synthetase_keff = 65.

    def __init__(self, id, model, amino_acid, RNA, codon):
        ProcessData.__init__(self, id, model)
        model.tRNA_data.append(self)
        self.codon = codon
        self.amino_acid = amino_acid
        self.RNA = RNA
        self.modifications = defaultdict(int)


class TranslocationData(ProcessData):
    """

    """

    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        self.keff = 65.
        # enzyme_dict = {enzyme_id: {length_dependent: <True or False>,
        #                            fixed_keff: <True or False>}}
        self.enzyme_dict = {}
        self.length_dependent_energy = False
        self.stoichiometry = {}
        model.translocation_data.append(self)


class PostTranslationData(ProcessData):
    """
    PostTranslationData id can be anything, but the preprocessed protein id
    and processed protein id must be defined
    """

    def __init__(self, id, model, processed_protein, preprocessed_protein):
        ProcessData.__init__(self, id, model)
        self.processed_protein_id = processed_protein
        self.unprocessed_protein_id = preprocessed_protein
        self.translocation = defaultdict(float)
        self.translocation_multipliers = {}

        self.folding_mechanism = ''
        self.aggregation_propensity = 0.
        # Dictionaries of {str(temperature): value}
        self.keq_folding = {}
        self.k_folding = {}
        # some small peptides are more likely to be folded by certain
        # chaperones. This is accounted for using propensity_scaling.
        self.propensity_scaling = 1.

        self.modifications = defaultdict(float)
        self.subreactions = defaultdict(float)
        self.surface_area = {}
        model.posttranslation_data.append(self)
