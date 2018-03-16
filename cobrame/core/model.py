from __future__ import print_function, division, absolute_import

import re
from six import iteritems
from warnings import warn

from cobra import Model, DictList
import numpy as np
from scipy.sparse import dok_matrix

from cobrame.core.reaction import (SummaryVariable, MetabolicReaction,
                                   TranscriptionReaction, TranslationReaction)
from cobrame.core.component import (Constraint, ProcessedProtein, Complex,
                                    TranslatedGene, TranscribedGene)
from cobrame.core import processdata
from cobrame.util import mu


class MEModel(Model):
    def __init__(self, *args):
        Model.__init__(self, *args)
        self.global_info = {}
        self.process_data = DictList()
        # create the biomass/dilution constraint
        self._biomass = Constraint("biomass")
        self._biomass_dilution = SummaryVariable("biomass_dilution")
        self._biomass_dilution.add_metabolites({self._biomass: -1})
        self.add_reactions([self._biomass_dilution])
        self._biomass_dilution.upper_bound = mu
        self._biomass_dilution.lower_bound = mu
        # maintenance energy
        self._gam = 0.
        self._ngam = 0.
        # Unmodeled protein is handled by converting protein_biomass to
        # biomass, and requiring production of the appropriate amount of dummy
        # protein
        self._unmodeled_protein_fraction = None

    def add_biomass_constraints_to_model(self, biomass_types):
        for biomass_type in biomass_types:
            if '_biomass' not in biomass_type:
                raise ValueError('Biomass types should be suffixed with '
                                 '"_biomass"')
            constraint_obj = Constraint(biomass_type)
            summary_variable_obj = SummaryVariable("%s_to_biomass" %
                                                   biomass_type)
            summary_variable_obj.add_metabolites({constraint_obj: -1,
                                                  self._biomass: 1})
            self.add_reactions([summary_variable_obj])

    @property
    def unmodeled_protein(self):
        return self.metabolites.get_by_id("protein_dummy")

    @property
    def unmodeled_protein_fraction(self):
        return self._unmodeled_protein_fraction

    @property
    def unmodeled_protein_biomass(self):
        return self.metabolites.get_by_id('unmodeled_protein_biomass')

    @unmodeled_protein_fraction.setter
    def unmodeled_protein_fraction(self, value):
        if 'protein_biomass_to_biomass' not in self.reactions:
            raise UserWarning("Must add SummaryVariable handling the protein"
                              "biomass constraint (via "
                              ":meth:`add_biomass_constraints_to_model`) "
                              "before defining the unmodeled protein fraction")

        # see the Biomass_formulations for an explanation
        amount = value / (1 - value)
        self.reactions.protein_biomass_to_biomass.add_metabolites(
            {self.unmodeled_protein_biomass: -amount}, combine=False)
        self.reactions.protein_biomass_to_biomass.add_metabolites(
            {self._biomass: 1 + amount}, combine=False)
        self._unmodeled_protein_fraction = value

    @property
    def gam(self):
        return self._gam

    @gam.setter
    def gam(self, value):
        if 'GAM' not in self.reactions:
            warn('Adding GAM reaction to model')
            self.add_reactions([SummaryVariable("GAM")])
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
            self.add_reactions([SummaryVariable("ATPM")])
            self.reactions.ATPM.add_metabolites(atp_hydrolysis)
        self.reactions.ATPM.lower_bound = value
        self._ngam = value

    @property
    def stoichiometric_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.StoichiometricData):
                yield data

    @property
    def complex_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.ComplexData):
                yield data

    @property
    def translation_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.TranslationData):
                yield data

    @property
    def transcription_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.TranscriptionData):
                yield data

    @property
    def generic_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.GenericData):
                yield data

    @property
    def tRNA_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.tRNAData):
                yield data

    @property
    def translocation_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.TranslocationData):
                yield data

    @property
    def posttranslation_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.PostTranslationData):
                yield data

    @property
    def subreaction_data(self):
        for data in self.process_data:
            if isinstance(data, processdata.SubreactionData):
                yield data

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

    def construct_s_matrix(self, growth_rate):
        """build the stoichiometric matrix at a specific growth rate"""
        # intialize to 0
        s = dok_matrix((len(self.metabolites), len(self.reactions)))
        # populate with stoichiometry
        for i, r in enumerate(self.reactions):
            for met, value in iteritems(r._metabolites):
                met_index = self.metabolites.index(met)
                if hasattr(value, "subs"):
                    s[met_index, i] = float(value.subs(mu, growth_rate))
                else:
                    s[met_index, i] = float(value)
        return s

    def construct_attribute_vector(self, attr_name, growth_rate):
        """build a vector of a reaction attribute at a specific growth rate

        Mainly used for upper and lower bounds"""
        return np.array([float(value.subs(mu, growth_rate))
                         if hasattr(value, "subs") else float(value)
                         for value in self.reactions.list_attr(attr_name)])

    def compute_solution_error(self, solution=None):
        errors = {}
        if solution is None:
            solution = self.solution
        s = self.construct_s_matrix(solution.f)
        lb = self.construct_attribute_vector("lower_bound", solution.f)
        ub = self.construct_attribute_vector("upper_bound", solution.f)
        x = np.array(solution.x)
        err = abs(s * x)
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

    def prune(self, skip=None):
        """remove all unused metabolites and reactions

        This should be run after the model is fully built. It will be
        difficult to add new content to the model once this has been run.

        skip: list
            List of complexes/proteins/mRNAs/TUs to remain unpruned from model.
        """
        if not skip:
            skip = []

        complex_data_list = [i.id for i in self.complex_data
                             if i.id not in skip]
        for c_d in complex_data_list:
            c = self.process_data.get_by_id(c_d)
            cplx = c.complex
            if len(cplx.reactions) == 1:
                list(cplx.reactions)[0].delete(remove_orphans=True)
                self.process_data.remove(self.process_data.get_by_id(c_d))

        for p in self.metabolites.query('_folded'):
            if 'partially' not in p.id and p.id not in skip:
                delete = True
                for rxn in p.reactions:
                    if rxn.metabolites[p] < 0:
                        delete = False
                        break

                if delete:
                    while len(p.reactions) > 0:
                        list(p.reactions)[0].delete(remove_orphans=True)
                        for data in self.process_data.query(p.id):
                            self.process_data.remove(data.id)

        for p in self.metabolites.query(re.compile('^protein_')):
            if isinstance(p, ProcessedProtein) and p.id not in skip:
                delete = True
                for rxn in p.reactions:
                    if rxn.metabolites[p] < 0:
                        delete = False
                        break
                if delete:
                    for rxn in list(p.reactions):
                        self.process_data.remove(rxn.posttranslation_data.id)
                        rxn.delete(remove_orphans=True)

        for p in self.metabolites.query(re.compile('^protein_')):
            if isinstance(p, TranslatedGene) and p.id not in skip:
                delete = True
                for rxn in p.reactions:

                    if rxn.metabolites[p] < 0 and not rxn.id.startswith(
                            'degradation'):
                        delete = False
                        break

                if delete:
                    for rxn in list(p.reactions):
                        p_id = p.id.replace('protein_', '')
                        data = self.process_data.get_by_id(p_id)
                        self.process_data.remove(data.id)
                        rxn.delete(remove_orphans=True)

        removed_rna = set()
        for m in list(self.metabolites.query(re.compile("^RNA_"))):

            delete = False if m.id in skip else True

            for rxn in m.reactions:
                if rxn.metabolites[m] < 0 and not rxn.id.startswith('DM_'):
                    delete = False
            if delete:
                try:
                    self.reactions.get_by_id('DM_' + m.id).remove_from_model(
                        remove_orphans=True)
                    if m in self.metabolites:
                        # Defaults to subtractive when removing reaction
                        m.remove_from_model()
                except KeyError:
                    pass
                else:
                    removed_rna.add(m.id)

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
                self.process_data.remove(t_process_id)
            else:
                # gets rid of the removed RNA from the products
                self.process_data.get_by_id(
                    t_process_id).RNA_products.difference_update(removed_rna)

            # update to update the TranscriptionReaction mRNA biomass
            # stoichiometry with new RNA_products
            if not delete:
                t.update()

    def remove_genes_from_model(self, gene_list):
        for gene in gene_list:
            # defaults to subtractive when removing model
            self.metabolites.get_by_id('RNA_'+gene).remove_from_model()
            protein = self.metabolites.get_by_id('protein_'+gene)
            for cplx in protein.complexes:
                print('Complex (%s) removed from model' % cplx.id)
                for rxn in cplx.metabolic_reactions:
                    try:
                        self.process_data.remove(rxn.id.split('_')[0])
                    except ValueError:
                        pass
                    rxn.remove_from_model()

            # If cannot import SymbolicParameter, assume using cobrapy
            # versions <= 0.5.11
            try:
                from optlang.interface import SymbolicParameter
            except ImportError:
                protein.remove_from_model(method='destructive')
            else:
                protein.remove_from_model(destructive=True)

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
                self.process_data.remove(t_process_id)

    def set_sasa_keffs(self, median_keff):
        # Get median SASA value considering all complexes in model
        sasa_list = []
        for met in self.metabolites:
            cplx_sasa = 0.
            if not isinstance(met, Complex):
                continue
            cplx_sasa += met.formula_weight ** (3. / 4)
            sasa_list.append(cplx_sasa)
        median_sasa = np.median(np.array(sasa_list))

        # redo scaling average SASA to 65.
        for rxn in self.reactions:
            if hasattr(rxn, 'keff') and rxn.complex_data is not None:
                sasa = rxn.complex_data.complex.formula_weight ** (3. / 4.)
                if sasa == 0:
                    raise UserWarning('No SASA for %s' % rxn)
                rxn.keff = sasa * median_keff / median_sasa
        for data in self.process_data:
            sasa = 0.
            if isinstance(data, processdata.TranslocationData):
                continue
            if hasattr(data, 'keff') and data.enzyme is not None:
                cplxs = \
                    [data.enzyme] if type(data.enzyme) == str else data.enzyme
                for cplx in cplxs:
                    sasa += \
                        self.metabolites.get_by_id(cplx).formula_weight ** \
                        (3. / 4)
                if sasa == 0:
                    raise UserWarning('No SASA for %s' % rxn)
                data.keff = sasa * median_keff / median_sasa

        self.update()
