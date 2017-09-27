from __future__ import print_function, division, absolute_import

from six import iteritems, string_types

from cobra import Metabolite as Component

from cobrame.util import dogma, mass


class MEComponent(Component):

    def __init__(self, id):
        Component.__init__(self, id)
        pass

    def remove_from_me_model(self, method='subtractive'):
        if self.id in self._model.process_data:
            self._model.process_data.remove(self.id)
        if self.id in self._model.complex_data:
            self._model.complex_data.remove(self.id)

        # If cannot import SymbolicParameter, assume using cobrapy
        # versions <= 0.5.11
        try:
            from optlang.interface import SymbolicParameter
        except ImportError:
            self.remove_from_model(method=method)
        else:
            if method.lower() == 'subtractive':
                self.remove_from_model(destructive=False)
            elif method.lower() == 'destructive':
                self.remove_from_model(destructive=True)
            else:
                raise AttributeError("method must be subtractive or "
                                     "destructive")


class Metabolite(MEComponent):
    pass


class TranscribedGene(MEComponent):

    def __init__(self, id):
        MEComponent.__init__(self, id)
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

    @property
    def mass(self):
        return mass.compute_rna_mass(self.nucleotide_sequence)


class TranslatedGene(MEComponent):
    @property
    def translation_data(self):
        locus = self.id.replace('protein_', '')
        return self._model.translation_data.get_by_id(locus)

    @property
    def complexes(self):
        """read-only link to the complexes that the gene forms"""
        complex_list = []
        for reaction in self.reactions:
            if reaction.__class__.__name__ == 'ComplexFormation':
                complex_list.append(reaction.complex)
        return complex_list

    @property
    def metabolic_reactions(self):
        """read-only link to the metabolic reactions that the gene is involved
        in"""
        metabolic_reactions = []
        for complexes in self.complexes:
            metabolic_reactions.extend(complexes.metabolic_reactions)
        return metabolic_reactions

    @property
    def mass(self):
        return self.translation_data.mass

    @property
    def amino_acid_sequence(self):
        return self.translation_data.amino_acid_sequence

    def get_surface_area(self, location):
        thickness = self._model.global_info['membrane_thickness'][location]
        mass = self.mass
        nm2_per_m2 = 1e18
        molecules_per_mmol = 6.022e20
        # return mass dependent SA from Liu et al 2014
        return 1.21 / thickness * 2. * mass * molecules_per_mmol / nm2_per_m2


class ProcessedProtein(MEComponent):

    @property
    def unprocessed_protein(self):
        return self._model.metabolites.get_by_id(self.unprocessed_protein_id)

    def __init__(self, id, unprocessed_protein_id):
        MEComponent.__init__(self, id)
        self.unprocessed_protein_id = unprocessed_protein_id

    @property
    def mass(self):
        return self.unprocessed_protein.mass


class Complex(MEComponent):
    @property
    def metabolic_reactions(self):
        """read-only link to MetabolicReactions"""
        reaction_list = []
        for reaction in self.reactions:
            if reaction.__class__.__name__ == 'MetabolicReaction':
                reaction_list.append(reaction)
        return reaction_list

    @property
    def mass(self):
        value = 0
        complex_data = self._model.process_data.get_by_id(self.id)
        for protein, coefficient in iteritems(complex_data.stoichiometry):
            protein_met = self._model.metabolites.get_by_id(protein)
            value += protein_met.mass * coefficient
        return value


class Ribosome(Complex):
    pass


class GenericComponent(MEComponent):
    pass


class GenerictRNA(MEComponent):
    pass


class RNAP(Complex):
    pass


class Constraint(MEComponent):
    pass


def create_component(component_id, default_type=MEComponent, rnap_set=set()):
    """creates a component and attempts to set the correct type"""
    if not isinstance(component_id, string_types):
        raise TypeError("%s must be a str, not %s" %
                        (repr(component_id), str(type(component_id))))
    if component_id.startswith("protein_"):
        return TranslatedGene(component_id)
    elif component_id.startswith("RNA_"):
        return TranscribedGene(component_id)
    elif component_id.startswith("ribosome"):
        return Ribosome(component_id)
    elif component_id.startswith("RNA_Polymerase") or component_id in rnap_set:
        return RNAP(component_id)
    elif component_id.startswith("generic_tRNA"):
        return GenerictRNA(component_id)
    elif component_id.endswith('_c'):
        return Component(component_id)
    else:
        return default_type(component_id)
