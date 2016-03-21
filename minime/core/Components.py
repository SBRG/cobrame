from cobra import Metabolite as Component
from minime.util import dogma
from six import iteritems


class Metabolite(Component):
    pass


class TranscribedGene(Component):

    def __init__(self, id):
        Component.__init__(self, id)
        self.left_pos = None
        self.right_pos = None
        self.strand = None
        self.RNA_type = ''
        self.nucleotide_sequence = ''

    @property
    def nucleotide_count(self):
        seq = self.nucleotide_sequence
        counts = {i: seq.count(i) for i in ("A", "T", "G", "C")}
        monophosphate_counts = {dogma.transcription_table[k].replace("tp_c",
                                                                     "mp_c"): v
                                for k, v in iteritems(counts)}
        return monophosphate_counts


class TranslatedGene(Component):
    @property
    def complexes(self):
        """read-only link to the complexes that the gene forms"""
        complex_list = []
        for reaction in self.reactions:
            if reaction.__class__.__name__ == 'ComplexFormation':
                complex_list.append(reaction.complex)
        return complex_list


class Complex(Component):
    @property
    def metabolic_reactions(self):
        """read-only link to MetabolicReactions"""
        reaction_list = []
        for reaction in self.reactions:
            if reaction.__class__.__name__ == 'MetabolicReaction':
                reaction_list.append(reaction)
        return reaction_list


class Ribosome(Complex):
    pass


class GenericComponent(Component):
    pass


class GenerictRNA(Component):
    pass


class RNAP(Complex):
    pass


class Constraint(Component):
    pass


class ModComplex(Component):
    pass


def create_component(component_id, default_type=Component):
    """creates a component and attempts to set the correct type"""
    if not isinstance(component_id, str):
        raise TypeError("%s must be a str, not %s" %
                        (repr(component_id), str(type(component_id))))
    if component_id.startswith("protein_"):
        return TranslatedGene(component_id)
    elif component_id.startswith("RNA_"):
        return TranscribedGene(component_id)
    elif component_id.startswith("ribosome"):
        return Ribosome(component_id)
    elif component_id.startswith("RNA_Polymerase"):
        return RNAP(component_id)
    elif component_id.startswith("generic_tRNA"):
        return GenerictRNA(component_id)
    else:
        return default_type(component_id)
