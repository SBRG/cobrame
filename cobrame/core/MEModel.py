from __future__ import print_function, division, absolute_import

import re

from cobra import Model, DictList
from numpy import array
from scipy.sparse import dok_matrix
from six import iteritems

from cobrame.core.Components import Constraint
from cobrame.core.ProcessData import *
from cobrame.core.MEReactions import *
from cobrame.util import mu


class MEModel(Model):
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
        # maintenance energy
        self._gam = 0.
        self._ngam = 0.
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
            self._DNA_biomass: -1,
            self._biomass: 1,
        })

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
        # see the Biomass_formulations for an explanation
        amount = value / (1 - value)
        self._protein_biomass_dilution.add_metabolites(
                {self.unmodeled_protein_biomass: -amount}, combine=False)
        self._protein_biomass_dilution.add_metabolites(
            {self._biomass: 1 + amount}, combine=False)
        self._unmodeled_protein_fraction = value

    @property
    def gam(self):
        return self._gam

    @gam.setter
    def gam(self, value):
        if 'GAM' not in self.reactions:
            warn('Adding GAM reaction to model')
            self.add_reaction(SummaryVariable("GAM"))
            self.reactions.GAM.lower_bound = mu
        atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'h_c': 1,
                          'pi_c': 1}
        for met, coeff in iteritems(atp_hydrolysis):
            self.reactions.GAM.add_metabolites({met: value * coeff},
                                               combine=False)
        self._gam = value

    @property
    def ngam(self):
        return self._ngam

    @ngam.setter
    def ngam(self, value):
        if 'ATPM' not in self.reactions:
            warn('Adding ATPM reaction to model')
            atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'h_c': 1,
                              'pi_c': 1}
            self.add_reaction(SummaryVariable("ATPM"))
            self.reactions.ATPM.add_metabolites(atp_hydrolysis)
        self.reactions.ATPM.lower_bound = value
        self._ngam = value

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
            for met, value in iteritems(r._metabolites):
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

    def prune(self,skip=[]):
        """remove all unused metabolites and reactions

        This should be run after the model is fully built. It will be
        difficult to add new content to the model once this has been run.

        skip: [str]
            List of complexes/proteins/mRNAs/TUs to remain unpruned from model.
        """
        complex_data_list = [i.id for i in self.complex_data if i.id not in skip]
        for c_d in complex_data_list:
            c = self.complex_data.get_by_id(c_d)
            cplx = c.complex
            if len(cplx.reactions) == 1:
                list(cplx.reactions)[0].delete(remove_orphans=True)
                self.complex_data.remove(self.complex_data.get_by_id(c_d))

        for p in self.metabolites.query('_folded'):
            if 'partially' not in p.id and p.id not in skip:
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
            if isinstance(p, ProcessedProtein) and p.id not in skip:
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
            if isinstance(p, TranslatedGene) and p.id not in skip:
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
            if m.id in skip:
                delete = False
            else:
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
            if t.id in skip:
                delete = False
            else:
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
                print('Complex (%s) removed from model' % cplx.id)
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

    def set_SASA_keffs(self, avg_keff):
        SASA_list = []
        for rxn in self.reactions:
            if hasattr(rxn, 'keff') and rxn.complex_data is not None:
                weight = rxn.complex_data.complex.mass
                SASA = weight ** (3. / 4.)
                if SASA == 0:
                    warn('Keff not updated for %s' % rxn)
                else:
                    SASA_list.append(SASA)

        for data in self.process_data:
            cplx_mw = 0.
            if isinstance(data, TranslocationData):
                continue
            if hasattr(data, 'keff') and data.enzyme is not None:
                if type(data.enzyme) == str:
                    cplxs = [data.enzyme]
                else:
                    cplxs = data.enzyme
                for cplx in cplxs:
                    if cplx in self.complex_data:
                        try:
                            cplx_mw += self.metabolites.get_by_id(cplx).mass \
                                       ** (3. / 4)
                        except:
                            warn('Complex (%s) cannot access mass' % cplx)
                    elif cplx.split('_mod_')[0] in self.complex_data:
                        cplx_mw += self.metabolites.get_by_id(cplx.split('_mod_')[0]).mass ** (3. / 4)

                SASA_list.append(cplx_mw)

        avg_SASA = array(SASA_list).mean()

        # redo scaling average SASA to 65.
        for rxn in self.reactions:
            if hasattr(rxn, 'keff') and rxn.complex_data is not None:
                weight = rxn.complex_data.complex.mass
                SASA = weight ** (3. / 4.)
                if SASA == 0:
                    SASA = avg_SASA
                rxn.keff = SASA * avg_keff / avg_SASA
        for data in self.process_data:
            SASA = 0.
            if isinstance(data, TranslocationData):
                continue
            if hasattr(data, 'keff') and data.enzyme is not None:
                if data.enzyme == str:
                    cplxs = [data.enzyme]
                else:
                    cplxs = data.enzyme
                for cplx in cplxs:
                    if cplx in self.complex_data:
                        SASA += self.metabolites.get_by_id(cplx).mass ** (3. / 4)
                    elif cplx.split('_mod_')[0] in self.complex_data:
                        SASA += self.metabolites.get_by_id(
                            cplx.split('_mod_')[0]).mass ** (3. / 4)

                if SASA == 0:
                    SASA = avg_SASA
                data.keff = SASA * avg_keff / avg_SASA

        self.update()
