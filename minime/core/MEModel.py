from numpy import array
from scipy.sparse import dok_matrix
import re

from cobra import Model, DictList

from minime.core.MEReactions import *
from minime.core.Components import Constraint
from minime.util import mu


class MEmodel(Model):
    def __init__(self, *args):
        Model.__init__(self, *args)
        self.global_info = {}
        self.stoichiometric_data = DictList()
        self.complex_data = DictList()
        self.modification_data = DictList()
        self.translation_data = DictList()
        self.transcription_data = DictList()
        self.generic_data = DictList()
        self.tRNA_data = DictList()
        self.translocation_data = DictList()
        self.posttranslation_data = DictList()
        self.subreaction_data = DictList()
        self.process_data = DictList()
        # create the biomass/dilution constraint
        self._biomass = Constraint("biomass")
        self._biomass_dilution = SummaryVariable("biomass_dilution")
        self._biomass_dilution.add_metabolites({self._biomass: -1})
        self.add_reaction(self._biomass_dilution)
        self._biomass_dilution.upper_bound = mu
        self._biomass_dilution.lower_bound = mu
        # Unmodeled protein is handled by converting protein_biomass to
        # biomass, and requiring production of the appropriate amount of dummy
        # protein
        self._unmodeled_protein_fraction = None
        self._protein_biomass = Constraint("protein_biomass")
        self._protein_biomass_dilution = SummaryVariable("protein_biomass_dilution")
        self._protein_biomass_dilution.add_metabolites({
            self._protein_biomass: -1,
            self._biomass: 1,
        })
        self._mRNA_biomass = Constraint("mRNA_biomass")
        self._mRNA_biomass_dilution = SummaryVariable("mRNA_biomass_dilution")
        self._mRNA_biomass_dilution.add_metabolites({
            self._mRNA_biomass: -1,
            self._biomass: 1,
        })
        self._tRNA_biomass = Constraint("tRNA_biomass")
        self._tRNA_biomass_dilution = SummaryVariable("tRNA_biomass_dilution")
        self._tRNA_biomass_dilution.add_metabolites({
            self._tRNA_biomass: -1,
            self._biomass: 1,
        })
        self._rRNA_biomass = Constraint("rRNA_biomass")
        self._rRNA_biomass_dilution = SummaryVariable("rRNA_biomass_dilution")
        self._rRNA_biomass_dilution.add_metabolites({
            self._rRNA_biomass: -1,
            self._biomass: 1,
        })

        self._ncRNA_biomass = Constraint("ncRNA_biomass")
        self._ncRNA_biomass_dilution = SummaryVariable("ncRNA_biomass_dilution")
        self._ncRNA_biomass_dilution.add_metabolites({
            self._ncRNA_biomass: -1,
            self._biomass: 1,
        })

        self._DNA_biomass = Constraint("DNA_biomass")
        self._DNA_biomass_dilution = SummaryVariable("DNA_biomass_dilution")
        self._DNA_biomass_dilution.add_metabolites({
            self._DNA_biomass: -1e-3,
            self._biomass: 1e-3,
        })

        self._DNA_biomass_dilution.lower_bound = mu
        self._DNA_biomass_dilution.upper_bound = mu

        self.add_reactions((self._protein_biomass_dilution,
                            self._mRNA_biomass_dilution,
                            self._tRNA_biomass_dilution,
                            self._rRNA_biomass_dilution,
                            self._ncRNA_biomass_dilution,
                            self._DNA_biomass_dilution))




    @property
    def unmodeled_protein(self):
        return self.metabolites.get_by_id("protein_dummy")

    @property
    def unmodeled_protein_fraction(self):
        return self._unmodeled_protein_fraction

    @property
    def unmodeled_protein_biomass(self):
        return self.metabolites.get_by_id('dummy_protein_biomass')

    @unmodeled_protein_fraction.setter
    def unmodeled_protein_fraction(self, value):
        # proportion = value / (1 - value)
        # see the Biomass_formulations for an explanation
        amount = value #/ self.unmodeled_protein.mass
        amount = value / (1 - value)
        self._protein_biomass_dilution.add_metabolites(
                {self.unmodeled_protein_biomass: -amount}, combine=False)
        self._protein_biomass_dilution.add_metabolites(
            {self._biomass: 1 + amount}, combine=False)
        self._unmodeled_protein_fraction = value

    def get_metabolic_flux(self, solution=None):
        """extract the flux state for metabolic reactions"""
        if solution is None:
            solution = self.solution
        if solution.status != "optimal":
            raise ValueError("solution status '%s' is not 'optimal'" %
                             solution.status)
        flux_dict = {r.id: 0 for r in self.stoichiometric_data}
        for reaction in self.reactions:
            if isinstance(reaction, MetabolicReaction):
                m_reaction_id = reaction.stoichiometric_data.id
                if reaction.reverse:
                    flux_dict[m_reaction_id] -= solution.x_dict[reaction.id]
                else:
                    flux_dict[m_reaction_id] += solution.x_dict[reaction.id]
            elif reaction.id.startswith("EX_") or reaction.id.startswith("DM"):
                flux_dict[reaction.id] = solution.x_dict[reaction.id]
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

    def compute_solution_error(self, solution=None):
        errors = {}
        if solution is None:
            solution = self.solution
        S = self.construct_S(solution.f)
        lb = self.construct_attribute_vector("lower_bound", solution.f)
        ub = self.construct_attribute_vector("upper_bound", solution.f)
        x = array(solution.x)
        err = abs(S * x)
        errors["max_error"] = err.max()
        errors["sum_error"] = err.sum()
        ub_err = min(ub - x)
        errors["upper_bound_error"] = abs(ub_err) if ub_err < 0 else 0
        lb_err = min(x - lb)
        errors["lower_bound_error"] = abs(lb_err) if lb_err < 0 else 0
        return errors

    def update(self):
        """updates all component reactions"""
        for r in self.reactions:
            if hasattr(r, "update"):
                r.update()

    def prune(self):
        """remove all unused metabolites and reactions

        This should be run after the model is fully built. It will be
        difficult to add new content to the model once this has been run.

        """
        complex_data_list = [i.id for i in self.complex_data]
        for c_d in complex_data_list:
            c = self.complex_data.get_by_id(c_d)
            cplx = c.complex
            if len(cplx.reactions) == 1:
                list(cplx.reactions)[0].delete(remove_orphans=True)
                self.complex_data.remove(self.complex_data.get_by_id(c_d))

        for p in self.metabolites.query('_folded'):
            if 'partially' not in p.id:
                delete = True
                for rxn in p._reaction:
                    try:
                        if p in rxn.reactants:
                            delete = False
                    except Exception as e:
                        print(rxn)
                        raise e
                if delete:
                    while len(p._reaction) > 0:
                        list(p._reaction)[0].delete(remove_orphans=True)
                        for data in self.posttranslation_data.query(p.id):
                            self.posttranslation_data.remove(data.id)

        for p in self.metabolites.query('partially_folded'):
            delete = True
            for rxn in p._reaction:
                try:
                    if p in rxn.reactants:
                        delete = False
                except Exception as e:
                    print(rxn)
                    raise e
            if delete:
                while len(p._reaction) > 0:
                    list(p._reaction)[0].delete(remove_orphans=True)
                    for data in self.posttranslation_data.query(p.id):
                        self.posttranslation_data.remove(data.id)

        for p in self.metabolites.query(re.compile('^protein_')):
            if isinstance(p, ProcessedProtein):
                delete = True
                for rxn in p._reaction:
                    try:
                        if p in rxn.reactants:
                            delete = False
                    except Exception as e:
                        print(rxn)
                        raise e
                if delete:
                    while len(p._reaction) > 0:
                        list(p._reaction)[0].delete(remove_orphans=True)
                        for data in self.posttranslation_data.query(p.id):
                            self.posttranslation_data.remove(data.id)

        for p in self.metabolites.query(re.compile('^protein_')):
            if isinstance(p, TranslatedGene):
                delete = True
                for rxn in p._reaction:
                    try:
                        if p in rxn.reactants and not rxn.id.startswith('degradation'):
                            delete = False
                    except Exception as e:
                        print(rxn)
                        raise e
                if delete:
                    while len(p._reaction) > 0:
                        list(p._reaction)[0].delete(remove_orphans=True)
                        p_id = p.id.replace('protein_', '')
                        for data in self.translation_data.query(p_id):
                            self.translation_data.remove(data.id)

        removed_RNA = set()
        for m in list(self.metabolites.query(re.compile("^RNA_"))):
            delete = True
            for rxn in m._reaction:
                if m in rxn.reactants and not rxn.id.startswith('DM_'):
                    delete = False
            if delete:
                try:
                    self.reactions.get_by_id('DM_' + m.id).remove_from_model(
                            remove_orphans=True)
                    if m in self.metabolites:
                        m.remove_from_model(method='subtractive')
                except KeyError:
                    pass
                else:
                    removed_RNA.add(m.id)

        for t in self.reactions.query('transcription_TU'):
            delete = True
            for product in t.products:
                if isinstance(product, TranscribedGene):
                    delete = False
            t_process_id = t.id.replace('transcription_', '')
            if delete:
                t.remove_from_model(remove_orphans=True)
                self.transcription_data.remove(t_process_id)
            else:
                # gets rid of the removed RNA from the products
                self.transcription_data.get_by_id(
                    t_process_id).RNA_products.difference_update(removed_RNA)

            # update to update the TranscriptionReaction mRNA biomass
            # stoichiometry with new RNA_products
            if not delete:
                t.update()

    def remove_genes_from_model(self, gene_list):
        for gene in gene_list:
            self.metabolites.get_by_id('RNA_'+gene).remove_from_model(method='subtractive')
            protein = self.metabolites.get_by_id('protein_'+gene)
            for cplx in protein.complexes:
                print cplx
                for rxn in cplx.metabolic_reactions:
                    try:
                        self.stoichiometric_data.remove(rxn.id.split('_')[0])
                    except ValueError:
                        pass
                    rxn.remove_from_model()

            protein.remove_from_model(method='destructive')

        # Remove all transcription reactions that now do not form a used
        # transcript
        for t in self.reactions.query('transcription_TU'):
            delete = True
            for product in t.products:
                if isinstance(product, TranscribedGene):
                    delete = False
            if delete:
                t.remove_from_model(remove_orphans=True)
                t_process_id = t.id.replace('transcription_', '')
                self.transcription_data.remove(t_process_id)

    def get_biomass_composition(self, solution=None):
        if solution is None:
            solution = self.solution
        biomass_composition = defaultdict(float)
        for met, stoich in self._protein_biomass_dilution.metabolites.items():
            if abs(stoich) > 1:
                protein_stoich = stoich
        biomass_composition['Protein'] = \
            solution.x_dict['protein_biomass_dilution'] * protein_stoich
        biomass_composition['tRNA'] = \
            solution.x_dict['tRNA_biomass_dilution']
        biomass_composition['mRNA'] = \
            solution.x_dict['mRNA_biomass_dilution']
        biomass_composition['ncRNA'] = \
            solution.x_dict['ncRNA_biomass_dilution']
        biomass_composition['rRNA'] = \
            solution.x_dict['rRNA_biomass_dilution']
        biomass_composition['Other'] = \
            solution.x_dict['biomass_component_dilution']

        return biomass_composition

    def RNA_to_protein_ratio(self, solution=None):
        composition = self.get_biomass_composition(solution=solution)
        RNA_to_protein = (composition['mRNA'] + composition['tRNA'] +
                          composition['rRNA'] + composition['ncRNA']) / \
                         (composition['Protein'] +
                          composition['Unmodeled Protein'])
        return RNA_to_protein

    def get_RNA_fractions_dict(self, solution=None):
        RNA_fractions = {}
        composition = self.get_biomass_composition(solution=solution)

        tRNA_to_RNA = (composition['tRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
        RNA_fractions['tRNA'] = tRNA_to_RNA

        rRNA_to_RNA = (composition['rRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
        RNA_fractions['rRNA'] = rRNA_to_RNA

        mRNA_to_RNA = (composition['mRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
        RNA_fractions['mRNA'] = mRNA_to_RNA

        ncRNA_to_RNA = (composition['ncRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
        RNA_fractions['ncRNA'] = ncRNA_to_RNA

        return RNA_fractions

    def get_membrane_composition(self, membrane='Inner_Membrane',
                                 solution=None):
        if not solution:
            solution = self.solution

        composition = {}

        composition['Dummy_protein'] = solution.x_dict[
                                           'SA_demand_dummy_protein_' + membrane]
        composition['Modeled_protein'] = solution.x_dict[
                                             'SA_demand_protein_' + membrane]
        composition['Lipid'] = solution.x_dict[
                                   'SA_demand_lipid_' + membrane]

        if membrane == 'Outer_Membrane':
            composition['LPS'] = solution.x_dict[
                                     'SA_demand_lps_' + membrane]
            composition['Lipoprotein'] = solution.x_dict[
                                             'SA_demand_lipoprotein_' + membrane]
        return composition

    def get_membrane_protein(self, membrane='Inner_Membrane', solution=None):
        if not solution:
            solution = self.solution

        composition = defaultdict(float)

        membrane_SA = self.metabolites.get_by_id('SA_protein_' + membrane)
        for rxn in membrane_SA.reactions:
            if membrane_SA in rxn.products:
                coeff = rxn._metabolites[membrane_SA]
                flux = coeff * solution.x_dict[rxn.id]
                if flux != 0:
                    composition[
                        rxn.id.replace('translocation_', '')] += flux
            elif 'lipid_modification' in rxn.id:
                coeff = rxn._metabolites[membrane_SA]
                flux = coeff * solution.x_dict[rxn.id]
                if flux != 0:
                    composition[rxn.id.split('_')[0]] += flux
        composition['dummy'] = solution.x_dict[
                                   'SA_demand_dummy_protein_' + membrane]
        return composition

    def make_composition_piechart(self, kind='Biomass', solution=None):
        try:
            import brewer2mpl
        except ImportError:
            color_map = None
        else:
            color_map = brewer2mpl.wesanderson.Zissou.mpl_colormap

        try:
            import pandas
        except ImportError:
            raise Exception("Pandas must be installed to get biomass piechart")

        if solution is None:
            solution = self.solution

        summary = {}
        if kind == 'Biomass':
            summary['Biomass composition'] = \
                self.get_biomass_composition(solution=solution)
            frame = pandas.DataFrame.from_dict(summary) / solution.f

        elif kind == 'Inner_Membrane':
            SA = solution.x_dict['SA_demand_Inner_Membrane']
            summary['Inner_Membrane'] = \
                self.get_membrane_composition(solution=solution, membrane=kind)
            frame = pandas.DataFrame.from_dict(summary) / SA
        elif kind == 'Outer_Membrane':
            SA = solution.x_dict['SA_demand_Outer_Membrane']
            summary['Outer_Membrane'] = \
                self.get_membrane_composition(solution=solution, membrane=kind)
            frame = pandas.DataFrame.from_dict(summary) / SA
        elif kind == 'Protein_Inner_Membrane':
            SA = solution.x_dict['SA_demand_protein_Inner_Membrane'] + \
                 solution.x_dict['SA_demand_dummy_protein_Inner_Membrane']
            summary['Inner_Membrane_Protein'] = \
                self.get_membrane_protein(solution=solution,
                                          membrane='Inner_Membrane')
            frame = pandas.DataFrame.from_dict(summary) / SA
        elif kind == 'Protein_Outer_Membrane':
            SA = solution.x_dict['SA_demand_protein_Outer_Membrane'] + \
                 solution.x_dict['SA_demand_dummy_protein_Outer_Membrane']
            summary['Outer_Membrane_Protein'] = \
                self.get_membrane_protein(solution=solution,
                                          membrane='Outer_Membrane')
            frame = pandas.DataFrame.from_dict(summary) / SA

        else:
            raise('%s not a valid composition kind' % kind)

        print('Component sum =', frame.sum().values[0])
        frame = frame[frame > 1e-4].dropna()
        return frame.plot(kind='pie', subplots=True, legend=None,
                          colormap=color_map)
