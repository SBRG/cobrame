from __future__ import print_function, absolute_import, division

from warnings import warn

from six import iteritems, string_types
from numpy import bool_, float_
from sympy import Basic, sympify, Symbol

import cobrame

_REQUIRED_REACTION_ATTRIBUTES = {"id", "name", "metabolites", "lower_bound",
                                 "upper_bound", "objective_coefficient",
                                 "variable_kind"}

# Reaction types can have different attributes
_REACTION_TYPE_DEPENDENCIES = \
    {'MetabolicReaction': ['complex_data',
                           'stoichiometric_data',
                           'keff', 'reverse'],
     'ComplexFormation': ['_complex_id',
                          'complex_data_id'],
     'PostTranslationReaction':
         ['posttranslation_data'],
     'TranscriptionReaction': ['transcription_data'],
     'GenericFormationReaction': [],
     'MEReaction': [],
     'SummaryVariable': [],
     'TranslationReaction': ['translation_data'],
     'tRNAChargingReaction': ['tRNA_data']}

_REQUIRED_PROCESS_DATA_ATTRIBUTES = {"id"}

# Process data types have different attributes
_PROCESS_DATA_TYPE_DEPENDENCIES = \
    {'StoichiometricData': ['_stoichiometry', 'lower_bound', 'upper_bound',
                            'subreactions'],

     'ComplexData': ['stoichiometry', 'complex_id', 'subreactions'],

     'TranscriptionData': ['subreactions', 'nucleotide_sequence',
                           'RNA_products', 'RNA_polymerase'],

     'TranslationData': ['subreactions', 'nucleotide_sequence', 'mRNA',
                         'protein'],

     'tRNAData': ['subreactions', 'codon', 'RNA', 'amino_acid',
                  'synthetase', 'synthetase_keff'],

     'TranslocationData': ['enzyme_dict', 'stoichiometry', 'keff',
                           'length_dependent_energy'],

     'PostTranslationData': ['processed_protein_id', 'unprocessed_protein_id',
                             'propensity_scaling', 'aggregation_propensity',
                             'translocation', 'subreactions', 'surface_area',
                             'keq_folding', 'k_folding', 'biomass_type',
                             'translocation_multipliers'],

     'SubreactionData': ['stoichiometry', 'enzyme', 'keff',
                         'element_contribution'],

     'GenericData': ['component_list']
     }

_REQUIRED_METABOLITE_ATTRIBUTES = {"id", "name", "formula", "compartment"}

_OPTIONAL_METABOLITE_ATTRIBUTES = {"charge", "_bound", "_constraint_sense"}

# Some metabolite types require additional attributes
_METABOLITE_TYPE_DEPENDENCIES = \
    {'TranscribedGene': ['left_pos', 'right_pos', 'strand', 'RNA_type',
                         'nucleotide_sequence'],
     'ProcessedProtein': ['unprocessed_protein_id']
     }


mu_temp = Symbol('mu')


def get_sympy_expression(value):
    """
    Return sympy expression from json string using sympify

    mu is assumed to be positive but using sympify does not apply this
    assumption. The mu symbol produced from sympify is replaced with
    cobrame's mu value to ensure the expression can be used in the model.

    Parameters
    ----------
    value : str
        String representation of mu containing expression

    Returns
    -------
    sympy expression
        Numeric representation of string with cobrame's mu symbol substituted

    """

    expression_value = sympify(value)
    return expression_value.subs(mu_temp, cobrame.mu)


def get_numeric_from_string(string):
    """

    Parameters
    ----------
    string : str
        String representation of numeric expression

    Returns
    -------
    float or sympy expression
        Numeric representation of string

    """
    try:
        return float(string)
    except ValueError:
        return get_sympy_expression(string)


def _fix_type(value):
    """convert possible types to str, float, and bool"""
    # Because numpy floats can not be pickled to json
    if isinstance(value, string_types):
        return str(value)
    if isinstance(value, float_):
        return float(value)
    if isinstance(value, bool_):
        return bool(value)
    if isinstance(value, set):
        return list(value)
    if isinstance(value, Basic):
        return str(value)
    if hasattr(value, 'id'):
        return str(value.id)
    # if value is None:
    #     return ''
    return value


def _reaction_to_dict(reaction):
    new_reaction = {key: _fix_type(getattr(reaction, key))
                    for key in _REQUIRED_REACTION_ATTRIBUTES
                    if key != 'metabolites'}

    reaction_type = reaction.__class__.__name__
    new_reaction['reaction_type'] = {}
    new_reaction['reaction_type'][reaction_type] = {}

    for attribute in _REACTION_TYPE_DEPENDENCIES.get(reaction_type, []):
        reaction_attribute = getattr(reaction, attribute)

        new_reaction['reaction_type'][reaction_type][attribute] = \
            _fix_type(reaction_attribute)

    # Add metabolites
    new_reaction['metabolites'] = {}
    for met, value in reaction.metabolites.items():
        new_reaction['metabolites'][met.id] = _fix_type(value)

    return new_reaction


def _process_data_to_dict(data):
    process_data_type = data.__class__.__name__

    new_data = {key: _fix_type(getattr(data, key))
                for key in _REQUIRED_PROCESS_DATA_ATTRIBUTES}

    new_data['process_data_type'] = {}
    new_data['process_data_type'][process_data_type] = {}
    new_process_data_type_dict = \
        new_data['process_data_type'][process_data_type]

    special_list = ['subreactions', 'stoichiometry', 'enzyme_dict',
                    'surface_area', 'keq_folding' 'k_folding']

    for attribute in _PROCESS_DATA_TYPE_DEPENDENCIES[process_data_type]:
        if attribute not in special_list:
            data_attribute = getattr(data, attribute)

            new_process_data_type_dict[attribute] = _fix_type(data_attribute)

        elif attribute == 'enzyme_dict':
            new_process_data_type_dict[attribute] = {}
            for cplx, values in getattr(data, attribute).items():
                new_process_data_type_dict[attribute][cplx] = {}
                for property, value in values.items():
                    new_process_data_type_dict[attribute][cplx][property] = \
                        _fix_type(value)
        else:
            new_process_data_type_dict[attribute] = {}
            for metabolite, coefficient in getattr(data, attribute).items():
                new_process_data_type_dict[attribute][metabolite] = \
                    _fix_type(coefficient)

    return new_data


def _metabolite_to_dict(metabolite):

    metabolite_type = metabolite.__class__.__name__

    new_metabolite = {key: _fix_type(getattr(metabolite, key))
                      for key in _REQUIRED_METABOLITE_ATTRIBUTES}

    # Some metabolites require additional information to construct working
    # ME-model
    new_metabolite['metabolite_type'] = {}
    new_metabolite['metabolite_type'][metabolite_type] = {}
    for attribute in _METABOLITE_TYPE_DEPENDENCIES.get(metabolite_type, []):
        metabolite_attribute = getattr(metabolite, attribute)
        new_metabolite['metabolite_type'][metabolite_type][attribute] = \
            metabolite_attribute

    return new_metabolite


def _get_attribute_array(dictlist, type):
    if type == 'reaction':
        return [_reaction_to_dict(reaction) for reaction in dictlist]
    elif type == 'process_data':
        return [_process_data_to_dict(data) for data in dictlist]
    elif type == 'metabolite':
        return [_metabolite_to_dict(metabolite) for metabolite in dictlist]
    else:
        raise TypeError('Type must be reaction, process_data or metabolite')


def _get_global_info_dict(global_info):
    new_global_info = {}
    for key, value in global_info.items():
        if type(value) != dict:
            new_global_info[key] = _fix_type(value)
        else:
            new_global_info[key] = value
    return new_global_info


def me_model_to_dict(model):
    """
    Create dictionary representation of full ME-model

    Parameters
    ----------
    model : :class:`~cobrame.core.model.MEModel`

    Returns
    -------
    dict
        Dictionary representation of ME-model

    """

    obj = dict(
        reactions=_get_attribute_array(model.reactions, 'reaction'),
        process_data=_get_attribute_array(model.process_data,
                                          'process_data'),
        metabolites=_get_attribute_array(model.metabolites, 'metabolite'),
        global_info=_get_global_info_dict(model.global_info)
    )

    return obj


# -----------------------------------------------------------------------------
# Functions below here are used to create a ME-model from its dictionary
# representation


def _add_metabolite_from_dict(model, metabolite_info):
    """
    Builds metabolite instances defined in dictionary, then add it to the
    ME-model being constructed.

    ProcessedProteins require additional information
    """

    metabolite_type_dict = metabolite_info['metabolite_type']
    if len(metabolite_type_dict) != 1:
        raise Exception('Only 1 metabolite_type in valid json')

    metabolite_type = list(metabolite_type_dict.keys())[0]

    # ProcessedProtein types require their unprocessed protein id as well
    if metabolite_type == 'ProcessedProtein':
        unprocessed_id = \
            metabolite_type_dict['ProcessedProtein']['unprocessed_protein_id']

        metabolite_obj = \
            getattr(cobrame, metabolite_type)(metabolite_info['id'],
                                              unprocessed_id)

    elif metabolite_type == 'TranscribedGene':
        rna_type = metabolite_type_dict['TranscribedGene']['RNA_type']
        nucleotide_sequence = \
            metabolite_type_dict['TranscribedGene']['nucleotide_sequence']
        metabolite_obj = \
            getattr(cobrame, metabolite_type)(metabolite_info['id'],
                                              rna_type, nucleotide_sequence)
    else:
        metabolite_obj = \
            getattr(cobrame, metabolite_type)(metabolite_info['id'])

    for attribute in _REQUIRED_METABOLITE_ATTRIBUTES:
        setattr(metabolite_obj, attribute, metabolite_info[attribute])

    for attribute in _METABOLITE_TYPE_DEPENDENCIES.get(metabolite_type, []):
        value = metabolite_type_dict[metabolite_type][attribute]
        setattr(metabolite_obj, attribute, value)

    model.add_metabolites([metabolite_obj])


def _add_process_data_from_dict(model, process_data_dict):
    """
    Builds process_data instances defined in dictionary, then add it to the
    ME-model being constructed.

    Most classes of process_data only require an id and model to initiate them,
    but TranslationData, tRNAData, PostTranslationData and GenericData require
    additional inputs.

    """

    # Create process data instances. Handel certain types individually
    id = process_data_dict['id']
    process_data_type_dict = process_data_dict['process_data_type']
    if len(process_data_type_dict) == 1:
        process_data_type, process_data_info = process_data_type_dict.popitem()
    else:
        print(process_data_type_dict, len(process_data_type_dict))
        raise Exception('Only 1 reaction_type in valid json')

    if process_data_type == 'TranslationData':
        mrna = process_data_info['mRNA']
        protein = process_data_info['protein']
        process_data = \
            getattr(cobrame, process_data_type)(id, model, mrna, protein)
    elif process_data_type == 'tRNAData':
        amino_acid = process_data_info['amino_acid']
        rna = process_data_info['RNA']
        codon = process_data_info['codon']
        process_data = \
            getattr(cobrame, process_data_type)(id, model, amino_acid, rna,
                                                codon)
    elif process_data_type == 'PostTranslationData':
        processed_protein_id = process_data_info['processed_protein_id']
        unprocessed_protein_id = process_data_info['unprocessed_protein_id']
        process_data = \
            getattr(cobrame, process_data_type)(id, model,
                                                processed_protein_id,
                                                unprocessed_protein_id)
    elif process_data_type == 'GenericData':
        component_list = process_data_info['component_list']
        process_data = \
            getattr(cobrame, process_data_type)(id, model, component_list)
        # Create reaction from generic process data
        process_data.create_reactions()
    else:
        process_data = getattr(cobrame, process_data_type)(id, model)

    # Set all of the required attributes using information in info dictionary
    for attribute in _REQUIRED_PROCESS_DATA_ATTRIBUTES:
        setattr(process_data, attribute, process_data_dict[attribute])

    # Some attributes depend on process data type. Set those here.
    for attribute in _PROCESS_DATA_TYPE_DEPENDENCIES.get(process_data_type,
                                                         []):
        value = process_data_info[attribute]
        try:
            setattr(process_data, attribute, value)
        except AttributeError:
            # set to the hidden attribute instead
            setattr(process_data, '_' + attribute, value)


def _add_reaction_from_dict(model, reaction_info):
    """
    Builds reaction instances defined in dictionary, then add it to the
    ME-model being constructed.

    """
    reaction_type_dict = reaction_info['reaction_type']

    if len(reaction_type_dict) == 1:
        reaction_type = list(reaction_type_dict.keys())[0]
        reaction_obj = getattr(cobrame, reaction_type)(reaction_info['id'])
    else:
        raise Exception('Only 1 reaction_type in valid json')

    for attribute in _REQUIRED_REACTION_ATTRIBUTES:
        # Metabolites are added to reactions using their update function,
        # skip setting metabolite stoichiometries here
        if attribute == 'metabolites':
            continue

        # upper and lower bounds may contain mu values. Handle that here
        value = reaction_info[attribute]
        if attribute in ['upper_bound', 'lower_bound']:
            value = get_sympy_expression(value)
        setattr(reaction_obj, attribute, value)

    # Some reactions are added to model when ME-models are initialized
    try:
        model.add_reactions([reaction_obj])
    except Exception:
        reaction_obj = model.reactions.get_by_id(reaction_obj.id)
        if reaction_type not in ['SummaryVariable',
                                 'GenericFormationReaction'] and \
                not reaction_obj.id.startswith('DM_'):
            warn('Reaction (%s) already in model' % reaction_obj.id)

    # These reactions types do not have update functions and need their
    # stoichiometries set explicitly .
    if reaction_type in ['SummaryVariable', 'MEReaction']:
        for key, value in reaction_info['metabolites'].items():
            reaction_obj.add_metabolites({key: get_sympy_expression(value)},
                                         combine=False)

    for attribute in _REACTION_TYPE_DEPENDENCIES.get(reaction_type, []):
        # Spontaneous reactions do no require complex_data
        if attribute == 'complex_data' and 'SPONT' in reaction_obj.id:
            continue

        value = reaction_type_dict[reaction_type][attribute]
        setattr(reaction_obj, attribute, value)

    if hasattr(reaction_obj, 'update'):
        reaction_obj.update()


def me_model_from_dict(obj):
    """
    Load ME-model from its dictionary representation. This will return
    a full :class:`~cobrame.core.model.MEModel` object identical to the
    one saved.

    Parameters
    ----------
    obj : dict
        Dictionary representation of ME-model

    Returns
    -------
    :class:`~cobrame.core.model.MEModel`:
        Full COBRAme ME-model
    """

    model = cobrame.MEModel()

    for k, v in iteritems(obj):
        if k in {'id', 'name', 'global_info'}:
            setattr(model, k, v)

    for metabolite in obj['metabolites']:
        _add_metabolite_from_dict(model, metabolite)

    for process_data in obj['process_data']:
        _add_process_data_from_dict(model, process_data)

    for reaction in obj['reactions']:
        _add_reaction_from_dict(model, reaction)

    model.update()

    return model
