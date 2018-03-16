from __future__ import print_function, division, absolute_import

from six import iteritems, string_types

from cobra import Metabolite as Component
from cobrame.util import dogma


class MEComponent(Component):
    """
    COBRAme component representation. Inherits from
    :class:`cobra.core.metabolite.Metabolite`

    Parameters
    ----------
    id : str
        Identifier of the component. Should follow best practices of child
        classes

    """
    def __init__(self, id):
        Component.__init__(self, id)

    def remove_from_me_model(self, method='subtractive'):
        """
        Remove metabolite from me model along with any relevant
        :class:`cobrame.core.processdata.ProcessData`

        Parameters
        ----------
        method : str
            - destructive: remove metabolite from model and remove reactions
              it is involved in
            - subtractive: remove only metabolite from model

        """
        if self.id in self._model.process_data:
            self._model.process_data.remove(self.id)

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
    """
    COBRAme metabolite representation

    Parameters
    ----------
    id : str
        Identifier of the metabolite

    """
    def __init__(self, id):
        MEComponent.__init__(self, id)


class TranscribedGene(MEComponent):
    """
    Metabolite class for gene created from
    :class:`cobrame.core.reaction.TranscriptionReaction`

    Parameters
    ----------
    id : str
        Identifier of the transcribed gene. As a best practice, this ID should
        be prefixed with 'RNA + _'

    RNA_type : str
        Type of RNA encoded by gene sequence (mRNA, rRNA, tRNA, or ncRNA)

    nucleotide_sequence : str
        String of base pair abbreviations for nucleotides contained in the gene

    Attributes
    ----------
    left_pos : int
        Left position of gene on the sequence of the (+) strain

    right_pos : int
        Right position of gene on the sequence of the (+) strain

    strand : str
        - (+) if the RNA product is on the leading strand
        - (-) if the RNA product is on the comple(mentary strand

    """

    def __init__(self, id, rna_type, nucleotide_sequence):
        MEComponent.__init__(self, id)
        self.left_pos = None
        self.right_pos = None
        self.strand = None
        self.RNA_type = rna_type
        self.nucleotide_sequence = nucleotide_sequence

    @property
    def nucleotide_count(self):
        """
        Get number of each nucleotide monophosphate

        Returns
        -------
        dict
            {nucleotide_monophosphate_id: count}

        """
        seq = self.nucleotide_sequence
        counts = {i: seq.count(i) for i in ("A", "T", "G", "C")}
        monophosphate_counts = {dogma.transcription_table[k].replace("tp_c",
                                                                     "mp_c"): v
                                for k, v in iteritems(counts)}
        return monophosphate_counts


class TranslatedGene(MEComponent):
    """
    Metabolite class for protein created from
    :class:`cobrame.core.reaction.TranslationReaction`

    Parameters
    ----------
    id : str
        Identifier of the translated protein product. Should be prefixed
        with "protein + _"

    """
    def __init__(self, id):
        MEComponent.__init__(self, id)

    @property
    def translation_data(self):
        """
        Get translation data that defines protein.

        Assumes that TranslatedGene is "protein + _ + <translation data id>"

        Returns
        -------
        :class:`cobrame.core.processdata.TranslationData`
            Translation data used to form translation reaction of protein
        """
        locus = self.id.replace('protein_', '')
        return self._model.process_data.get_by_id(locus)

    @property
    def complexes(self):
        """Get the complexes that the protein forms

        Returns
        -------
        list
            List of :class:`cobrame.core.component.Complex` s that the protein
            is a subunit of
        """
        complex_list = []
        for reaction in self.reactions:
            if hasattr(reaction, 'complex'):
                complex_list.append(reaction.complex)
        return complex_list

    @property
    def metabolic_reactions(self):
        """Get the mtabolic reactions that the protein helps catalyze

        Returns
        -------
        list
            List of :class:`cobrame.core.reactions.MetabolicReaction` s
            that the protein helps catalyze
        """
        metabolic_reactions = []
        for complexes in self.complexes:
            metabolic_reactions.extend(complexes.metabolic_reactions)
        return metabolic_reactions

    @property
    def amino_acid_sequence(self):
        """
        Get amino acid sequence of protein

        Returns
        -------
        str
            Amino acid sequence of protein

        """
        return self.translation_data.amino_acid_sequence


class ProcessedProtein(MEComponent):
    """
    Metabolite class for protein created from
    :class:`cobrame.core.reaction.PostTranslationReaction`

    Parameters
    ----------
    id : str
        Identifier of the processed protein

    unprocessed_protein_id : str
        Identifier of protein before being processed by PostTranslationReaction

    """
    def __init__(self, id, unprocessed_protein_id):
        MEComponent.__init__(self, id)
        self.unprocessed_protein_id = unprocessed_protein_id

    @property
    def unprocessed_protein(self):
        """Get unprocessed protein reactant in PostTranslationReaction

        Returns
        -------
        :class:`cobrame.core.component.TranslatedGene`
            Unprocessed protein object
        """
        return self._model.metabolites.get_by_id(self.unprocessed_protein_id)


class Complex(MEComponent):
    """
    Metabolite class for protein complexes

    Parameters
    ----------
    id : str
        Identifier of the protein complex.
    """

    def __init__(self, id):
        MEComponent.__init__(self, id)

    @property
    def metabolic_reactions(self):
        """Get metabolic reactions catalyzed by complex

        Returns
        -------
        list
            List of :class:`cobrame.core.reaction.MetabolicReaction` s
            catalyzed by complex.
        """
        reaction_list = []
        for reaction in self.reactions:
            if hasattr(reaction, 'stoichiometric_data'):
                reaction_list.append(reaction)
        return reaction_list


class Ribosome(Complex):
    """
    Metabolite class for Ribosome complexes. Inherits from
    :class:`cobrame.core.component.Complex`

    Parameters
    ----------
    id : str
        Identifier of the Ribosome.
    """
    def __init__(self, id):
        Complex.__init__(self, id)


class GenericComponent(MEComponent):
    """
    Metabolite class for generic components created from
    :class:`cobrame.core.reaction.GenericFormationReaction`

    Parameters
    ----------
    id : str
        Identifier of the generic tRNA. As a best practice should follow
        template: 'generic + _ +  <generic metabolite id>'
    """
    def __init__(self, id):
        MEComponent.__init__(self, id)


class GenerictRNA(MEComponent):
    """
    Metabolite class for generic tRNAs created from
    :class:`cobrame.core.reaction.tRNAChargingReaction`

    Parameters
    ----------
    id : str
        Identifier of the generic tRNA. As a best practice should follow
        template: 'generic_tRNA + _ + <codon> + _ + <amino acid metabolite id>'

    """
    def __init__(self, id):
        MEComponent.__init__(self, id)


class RNAP(Complex):
    """
    Metabolite class for RNA polymerase complexes. Inherits from
    :class:`cobrame.core.component.Complex`

    Parameters
    ----------
    id : str
        Identifier of the RNA Polymerase.
    """

    def __init__(self, id):
        Complex.__init__(self, id)


class Constraint(MEComponent):
    """
    Metabolite class for global constraints such as biomass

    Parameters
    ----------
    id : str
        Identifier of the constraint
    """
    def __init__(self, id):
        MEComponent.__init__(self, id)


def create_component(component_id, default_type=MEComponent, rnap_set=set()):
    """creates a component and attempts to set the correct type"""
    if not isinstance(component_id, string_types):
        raise TypeError("%s must be a str, not %s" %
                        (repr(component_id), str(type(component_id))))
    if component_id.startswith("protein_"):
        return TranslatedGene(component_id)
    elif component_id.startswith("RNA_"):
        raise ValueError('TranscribedGene (%s) should not be added using '
                         'create_component. It requires additional information'
                         'when creating instance.' % component_id)
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
