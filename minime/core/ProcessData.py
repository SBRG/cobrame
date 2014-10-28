from minime.util.dogma import *


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

    @property
    def model(self):
        return self._model

    @property
    def parent_reactions(self):
        return {self._model.reactions.get_by_id(i) for i in self._parent_reactions}

    def _update_parent_reactions(self):
        reactions = self._model.reactions
        for i in self._parent_reactions:
            reactions.get_by_id(i).update()

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))

class MetabolicReactionData(ProcessData):
    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.metabolic_reaction_data.append(self)

class ComplexData(ProcessData):
    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.complex_data.append(self)

class TranscriptionData(ProcessData):
    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.transcription_data.append(self)
        self.nucleotide_sequence=''


class TranslationData(ProcessData):
    protein_per_mRNA = 10

    def __init__(self, id, model):
        ProcessData.__init__(self, id, model)
        model.translation_data.append(self)

    def compute_sequence_from_DNA(self, dna_sequence):
        codons = (dna_sequence[i: i + 3] for i in range(0, (len(dna_sequence)), 3))
        self.amino_acid_sequence = ''.join(codon_table[i] for i in codons).rstrip("*")
