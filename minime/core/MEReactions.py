from cobra import Reaction
from collections import defaultdict

from minime.util.dogma import *
from minime.util import mu

class MetabolicReaction(Reaction):
    """Lumped metabolic reaction with enzyme complex formation"""

    # process data
    _complex_data = None
    @property
    def complex_data(self):
        return self._complex_data
    @complex_data.setter
    def complex_data(self, process_data):
        self._complex_data = process_data
        process_data._parent_reactions.add(self.id)

    _metabolic_reaction_data = None
    @property
    def metabolic_reaction_data(self):
        return self._metabolic_reaction_data
    @metabolic_reaction_data.setter
    def metabolic_reaction_data(self, process_data):
        self._metabolic_reaction_data = process_data
        process_data._parent_reactions.add(self.id)

    keff = 65.
    reverse = False

    def update(self):
        new_stoichiometry = {}
        for component_name, value in self.complex_data._stoichiometry.iteritems():
            component = self._model.metabolites.get_by_id(component_name)
            new_stoichiometry[component] = -value * mu / self.keff
        sign = -1 if self.reverse else 1
        for component_name, value in self.metabolic_reaction_data._stoichiometry.iteritems():
            component = self._model.metabolites.get_by_id(component_name)
            if component not in new_stoichiometry:
                new_stoichiometry[component] = 0
            new_stoichiometry[component] += value * sign
        # replace old stoichiometry with new one
        # TODO prune out old metabolites
        # doing checks for relationship every time is unnecessary. don't do.
        self.add_metabolites(new_stoichiometry, combine=False, add_to_container_model=False)


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
        new_stoichiometry = defaultdict(lambda: 0)
        try:
            transcript = self._model.metabolites.get_by_id("mRNA_" + protein_id)
            new_stoichiometry[transcript] = 1
        except KeyError:
            warn("TranscriptionRection 'mRNA_%s' not found in model" % protein_id)

        NT_mapping = {key: self._model.metabolites.get_by_id(value)
                      for key, value in transcription_table.iteritems()}
        for i in self.transcription_data.nucleotide_sequence:
            new_stoichiometry[NT_mapping[i]] -= 1
        self.add_metabolites(new_stoichiometry, combine=False, add_to_container_model=False)

class TranslationReaction(Reaction):
    """Translation of a TranscribedGene to a TranslatedGene"""
    _translation_data = None
    @property
    def translation_data(self):
        return self._translation_data
    @translation_data.setter
    def translation_data(self, process_data):
        self._translation_data = process_data
        process_data._parent_reactions.add(self.id)


    def update(self):
        protein_id = self.translation_data.id
        new_stoichiometry = defaultdict(lambda: 0)
        try:
            transcript = self._model.metabolites.get_by_id("mRNA_" + protein_id)
            new_stoichiometry[transcript] = -1./self.translation_data.protein_per_mRNA
        except KeyError:
            warn("TranscriptionRection 'mRNA_%s' not found in model" % protein_id)
        protein = self._model.metabolites.get_by_id("protein_" + protein_id)
        aa_mapping = {key: self._model.metabolites.get_by_id(value)
                      for key, value in amino_acids.iteritems()}
        for i in self.translation_data.amino_acid_sequence:
            new_stoichiometry[aa_mapping[i]] -= 1
        new_stoichiometry[protein] = 1
        self.add_metabolites(new_stoichiometry, combine=False, add_to_container_model=False)
