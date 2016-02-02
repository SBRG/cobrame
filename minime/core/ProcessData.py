from collections import defaultdict

from minime.util.dogma import *
from minime.util.mass import *
from minime.util import mu
from minime.core.MEReactions import *


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
    """Encodes the stoichiometry  for a reaction.

    Used by Metabolic Reactions
    """
    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.stoichiometric_data.append(self)
        self._stoichiometry = {}

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

    def get_complex_data(self):
        for i in self._model.complex_data:
            if self.id in i.modifications:
                yield i


class ComplexData(ProcessData):

    @property
    def formation(self):
        """a read-only link to the formation reaction"""
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
        return self.id if self._complex_id is None else self._complex_id

    @complex_id.setter
    def complex_id(self, value):
        self._complex_id = None if value == self.id else value

    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.complex_data.append(self)
        # {Component.id: stoichiometry}
        self.stoichiometry = defaultdict(float)
        self.translocation = {}
        self.chaperones = {}
        # {ModificationData.id : number}
        self.modifications = {}
        self._complex_id = None  # assumed to be the same as id if None

    def create_complex_formation(self, verbose=True):
        """creates a complex formation reaction

        This assumes none exists already. Will create a reaction (prefixed by
        'formation_') which forms the complex"""
        formation_id = "formation_" + self.id
        if formation_id in self._model.reactions:
            raise ValueError("reaction %s already in model" % formation_id)
        formation = ComplexFormation(formation_id)
        formation._complex_id = self.complex_id
        self._model.add_reaction(formation)
        formation.update(verbose=verbose)


class TranscriptionData(ProcessData):
    def __init__(self, id, model, RNA_products=set()):
        ProcessData.__init__(self, id, model)
        model.transcription_data.append(self)
        self.nucleotide_sequence = ''
        self.RNA_products = RNA_products
        # i.e. {"amp_c": 10, "gmp_c": 11, "ump_c": 9, "cmp_c": 11}
        self.excised_bases = {}
        # {ModificationData.id : number}
        self.modifications = defaultdict(int)
        # Used if not creating a "MiniME" model
        self.using_RNAP = True

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
    protein_per_mRNA = .5  # per second
    amino_acid_sequence = ""
    mRNA = None
    nucleotide_sequence = ""

    def __init__(self, id, model, mRNA, protein):
        ProcessData.__init__(self, id, model)
        model.translation_data.append(self)
        self.mRNA = mRNA
        self.protein = protein
        # Used if not creating a "MiniME" model
        self.using_ribosome = True
        self.subreactions = defaultdict(int)


    @property
    def amino_acid_sequence(self):
        if len(self.nucleotide_sequence) % 3 != 0:
            self.nucleotide_sequence = self.nucleotide_sequence[:-1]
            print self.id, ' Needs frameshift?'
        codons = (self.nucleotide_sequence[i: i + 3]
                  for i in range(0, (len(self.nucleotide_sequence)), 3))
        amino_acid_sequence = ''.join(codon_table[i] for i in codons)
        amino_acid_sequence = amino_acid_sequence.rstrip("*")
        if "*" in amino_acid_sequence:
            amino_acid_sequence = amino_acid_sequence.replace('*', 'U')
        return amino_acid_sequence

    @property
    def codon_count(self):
        # exclude the last three stop codons from count
        codons = (self.nucleotide_sequence[i: i + 3]
                  for i in range(0, (len(self.nucleotide_sequence)-3), 3))
        codon_count = defaultdict(int)
        for i in codons:
            if len(i) % 3 != 0:
                print self.id, 'Needs Frameshift?'
                continue
            codon_count[i.replace('T', 'U')] += 1

        # Remove one methionine (AUG) from codon count to account for start

        if self.nucleotide_sequence.startswith('ATG'):
            codon_count['AUG'] -= 1
        elif self.nucleotide_sequence.startswith('GTG'):
            codon_count['GUG'] -= 1
        elif self.nucleotide_sequence.startswith('TTG'):
            codon_count['UUG'] -= 1
        elif self.nucleotide_sequence.startswith('ATT'):
            codon_count['AUU'] -= 1
        elif self.nucleotide_sequence.startswith('CTG'):
            codon_count['CUG'] -= 1
        else:
            raise NameError('No start codon in DNA sequence %s' %self.mRNA)
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


class tRNAData(ProcessData):
    synthetase = None
    synthetase_keff = 65.

    def __init__(self, id, model, amino_acid, RNA, codon):
        ProcessData.__init__(self, id, model)
        model.tRNA_data.append(self)
        self.codon = codon
        self.amino_acid = amino_acid
        self.RNA = RNA
        self.tRNA_keff = 2.39 * mu / (mu + 0.391)


class ProteinTranslocationData(ProcessData):
    """
    The user will have to use update after all complexes are into the
    model...because otherwise how can you be sure that the translocases
    are added into the model?

    """

    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        self.keff = 0.
        self.costs_complexes = []
        model.translocation_pathways.append(self)

    def add_translocation_cost(self, model, complex, protein):
        """really need to document this

           Need to add in membrane area calculations
        """

        costs = {}

        def sec_translocation(protein):  # s
            """

            """

            protein, stoichiometry = protein.split('(')
            stoichiometry = float(stoichiometry[:-1])
            protein_length = len(
                model.translation_data.get_by_id(protein).amino_acid_sequence)
            atp = int(protein_length/25)

            # SecB
            costs[self.costs_complexes[0]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600
            # SecA
            costs[self.costs_complexes[1]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600
            # Sec-CPLX
            costs[self.costs_complexes[2]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600
            costs['atp_c'] = atp
            costs['adp_c'] = -atp
            costs['pi_c'] = -atp

        def tat_translocation(protein):  # t
            """

            """

            peptide, stoichiometry = protein.split('(')
            # e.g. (1:14) for 1 of this protein in the final complex:14
            # tatAs required
            s1, s2 = stoichiometry[:-1].split(':')

            # TatBC
            costs[self.costs_complexes[0]] = \
                mu * float(s1) / self.keff / 3600
            # TatA
            costs[self.costs_complexes[1]] = \
                mu * float(s1) * float(s2) / self.keff / 3600

        def bam_translocation(protein):  # b
            """

            """

            protein, stoichiometry = protein.split('(')
            stoichiometry = float(stoichiometry[:-1])

            # SurA
            costs[self.costs_complexes[0]] = \
                mu * stoichiometry / self.keff / 3600
            # Bam
            costs[self.costs_complexes[1]] = \
                mu * stoichiometry / self.keff / 3600

        def lol_translocation(protein):  # l
            """

            """

            protein, stoichiometry = protein.split('(')
            stoichiometry = float(stoichiometry[:-1])

            # LolCDE
            costs[self.costs_complexes[0]] = \
                mu * stoichiometry / self.keff / 3600
            # LolA
            costs[self.costs_complexes[1]] = \
                mu * stoichiometry / self.keff / 3600
            # LolB
            costs[self.costs_complexes[2]] = \
                mu * stoichiometry / self.keff / 3600

            costs['atp_c'] = 1
            costs['adp_c'] = -1
            costs['pi_c'] = -1

        def yidC_translocation(protein):  # y
            """

            """

            protein, stoichiometry = protein.split('(')
            stoichiometry = float(stoichiometry[:-1])
            protein_length = len(
                model.translation_data.get_by_id(protein).amino_acid_sequence)

            # SRP
            costs[self.costs_complexes[0]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600
            # YidC
            costs[self.costs_complexes[1]] = \
                mu * 2 / (self.keff/protein_length) / 3600

            costs['gtp_c'] = 1
            costs['gdp_c'] = -1
            costs['pi_c'] = -1

        def secA_translocation(protein):  # a
            """

            """

            # Not sure what the actualy number is...will need to do
            # calculations. For now, 1/3 of the peptide requires secA...since
            # secA binding isn't detected until after ribosome release
            # -Jo 10/31/13

            protein, stoichiometry = protein.split('(')
            stoichiometry = float(stoichiometry[:-1])
            protein_length = len(
                model.translation_data.get_by_id(protein).amino_acid_sequence)
            atp = int(protein_length / 3 / 25)

            # SecA
            costs[self.costs_complexes[0]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600

            costs['atp_c'] = atp
            costs['adp_c'] = -atp
            costs['pi_c'] = -atp

        def srp_yidC_translocation(protein):  # p
            """

            """

            protein, stoichiometry = protein.split('(')
            stoichiometry = float(stoichiometry[:-1])
            protein_length = len(
                model.translation_data.get_by_id(protein).amino_acid_sequence)

            # SRP
            costs[self.costs_complexes[0]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600
            # YidC
            costs[self.costs_complexes[1]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600
            # Sec
            costs[self.costs_complexes[2]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600

            costs['gtp_c'] = 1
            costs['gdp_c'] = -1
            costs['pi_c'] = -1

        def srp_translocation(protein):  # r
            """	FtsY keff needs to fixed"""
            protein, stoichiometry = protein.split('(')
            stoichiometry = float(stoichiometry[:-1])
            protein_length = len(
                model.translation_data.get_by_id(protein).amino_acid_sequence)

            # SRP
            costs[self.costs_complexes[0]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600
            # FtsY
            costs[self.costs_complexes[1]] = \
                mu * stoichiometry / 65. / 3600
            # Sec
            costs[self.costs_complexes[2]] = \
                mu * stoichiometry / (self.keff/protein_length) / 3600

            costs['gtp_c'] = 2
            costs['gdp_c'] = -2
            costs['pi_c'] = -2

        eval(self.id+'(protein)')
        formation = model.complex_data.get_by_id(complex).formation
        for metabolite, cost in costs.iteritems():
            model.complex_data.get_by_id(complex).\
                stoichiometry[metabolite] += cost
            # if membrane...then add into self.membrane
        formation.update()
