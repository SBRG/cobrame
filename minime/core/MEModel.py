from numpy import array
from scipy.sparse import dok_matrix

from cobra import Model, DictList

from minime.core.MEReactions import *
from minime.util import mu


class MEmodel(Model):
    def __init__(self, *args):
        Model.__init__(self, *args)
        self.metabolic_reaction_data = DictList()
        self.complex_data = DictList()
        self.translation_data = DictList()
        self.transcription_data = DictList()
        self.tRNA_data = DictList()
        # create the biomass/dilution constraint
        self._biomass = Constraint("biomass")
        self._biomass_dilution = SummaryVariable("biomass_dilution")
        self._biomass_dilution.add_metabolites({self._biomass: -1})
        self.add_reaction(self._biomass_dilution)
        self._biomass_dilution.upper_bound = mu
        self._biomass_dilution.lower_bound = mu

    def get_metabolic_flux(self, solution=None):
        """extract the flux state for metabolic reactions"""
        if solution is None:
            solution = self.solution
        if solution.status != "optimal":
            raise ValueError("solution status '%s' is not 'optimal'" %
                             solution.status)
        flux_dict = {r.id: 0 for r in self.metabolic_reaction_data}
        for reaction in self.reactions:
            if isinstance(reaction, MetabolicReaction):
                m_reaction_id = reaction.metabolic_reaction_data.id
                if reaction.reverse:
                    flux_dict[m_reaction_id] -= solution.x_dict[reaction.id]
                else:
                    flux_dict[m_reaction_id] += solution.x_dict[reaction.id]
        return flux_dict

    def get_transcription_flux(self, solution=None):
        """extract the transcription flux state"""
        if solution is None:
            solution = self.solution
        if solution.status != "optimal":
            raise ValueError("solution status '%s' is not 'optimal'" %
                             solution.status)
        flux_dict = {}
        for reaction in self.reactions:
            if isinstance(reaction, TranscriptionReaction):
                for rna_id in reaction.transcription_data.RNA_products:
                    locus_id = rna_id.replace("RNA_", "", 1)
                    if locus_id not in flux_dict:
                        flux_dict[locus_id] = 0
                    flux_dict[locus_id] += solution.x_dict[reaction.id]
        return flux_dict

    def get_translation_flux(self, solution=None):
        """extract the translation flux state"""
        if solution is None:
            solution = self.solution
        if solution.status != "optimal":
            raise ValueError("solution status '%s' is not 'optimal'" %
                             solution.status)
        flux_dict = {r.id: 0 for r in self.translation_data}
        for reaction in self.reactions:
            if isinstance(reaction, TranslationReaction):
                protein_id = reaction.translation_data.id
                flux_dict[protein_id] += solution.x_dict[reaction.id]
        return flux_dict

    def construct_S(self, growth_rate):
        """build the stoichiometric matrix at a specific growth rate"""
        # intialize to 0
        S = dok_matrix((len(self.metabolites), len(self.reactions)))
        # populate with stoichiometry
        for i, r in enumerate(self.reactions):
            for met, value in r._metabolites.iteritems():
                met_index = self.metabolites.index(met)
                if hasattr(value, "subs"):
                    S[met_index, i] = float(value.subs(mu, growth_rate))
                else:
                    S[met_index, i] = float(value)
        return S

    def construct_attribute_vector(self, attr_name, growth_rate):
        """build a vector of a reaction attribute at a specific growth rate

        Mainly used for upper and lower bounds"""
        return array([float(value.subs(mu, growth_rate))
                      if hasattr(value, "subs") else float(value)
                      for value in self.reactions.list_attr(attr_name)])
