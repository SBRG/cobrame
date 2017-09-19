import copy
import json
from warnings import warn
import os

from sympy import Basic, sympify, Symbol
from numpy import bool_, float_
from six import iteritems, string_types
from jsonschema import validate, ValidationError

from cobra.io.json import save_json_model, load_json_model
import cobrame
from cobrame import mu

mu_temp = Symbol('mu')

cur_dir = os.path.dirname(os.path.abspath(__file__))


def save_json_me(me0, file_name, pretty=False):
    """
    Save a stripped-down JSON version of the ME-model. This will exclude all of
    ME-Model information except the reaction stoichiometry information and the
    reaction bounds. Saving/loading a model in this format will thus occur much
    quicker, but limit the ability to edit the model and use most of its
    features.

    :param :class:`~cobrame.core.MEModel.MEModel` model:
        A full ME-model

    :param str or file-like object file_name:
        Filename of the JSON output

    :returns JSON-object:
        Stripped-down JSON representation of full ME-model
    """
    me = copy.deepcopy(me0)

    for rxn in me.reactions:
        for met in rxn.metabolites:
            s = rxn._metabolites[met]
            if isinstance(s, Basic):
                rxn._metabolites[met] = str(s)
        if isinstance(rxn.lower_bound, Basic):
            rxn.lower_bound = str(rxn.lower_bound)
        if isinstance(rxn.upper_bound, Basic):
            rxn.upper_bound = str(rxn.upper_bound)

    for met in me.metabolites:
        if isinstance(met._bound, Basic):
            met._bound = str(met._bound)

    save_json_model(me, file_name)


def get_sympy_expression(value):
    """
    Return sympy expression from json string using sympify


    mu is assumed to be positive but using sympify does not apply this
    assumption"""

    expression_value = sympify(value)
    return expression_value.subs(mu_temp, mu)


def load_json_me(file_name):
    """
    Load a stripped-down JSON version of the ME-model. This will exclude all of
    ME-Model information except the reaction stoichiometry information and the
    reaction bounds. Saving/loading a model in this format will thus occur much
    quicker, but limit the ability to edit the model and use most of its
    features.

    :param str or file-like object file_name:
        Filename of the JSON ME-model

    :returns COBRA Model object:
        COBRA Model representation of the full ME-model. This will not include
        all of the functionality of a :class:`~cobrame.core.MEModel.MEModel`
    """
    me = load_json_model(file_name)

    # Re-convert stoichiometries back to sympy
    for rxn in me.reactions:
        for met in rxn.metabolites:
            s = rxn._metabolites[met]
            try:
                rxn._metabolites[met] = float(s)
            except ValueError:
                rxn._metabolites[met] = get_sympy_expression(s)
        try:
            rxn.lower_bound = float(rxn.lower_bound)
        except ValueError:
            rxn.lower_bound = get_sympy_expression(rxn.lower_bound)
        try:
            rxn.upper_bound = float(rxn.upper_bound)
        except ValueError:
            rxn.upper_bound = get_sympy_expression(rxn.upper_bound)

    for met in me.metabolites:
        b = met._bound
        try:
            met._bound = float(b)
        except ValueError:
            met._bound = get_sympy_expression(b)

    return me

# -----------------------------------------------------------------------------
# Functions below here facilitate json dumping/loading of full ME-models with
# all process_data/reaction info intact.
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
     'Reaction': [],
     'SummaryVariable': [],
     'TranslationReaction': ['translation_data'],
     'tRNAChargingReaction': ['tRNA_data']}

_REQUIRED_PROCESS_DATA_ATTRIBUTES = {"id"}

# Process data types have different attributes
_PROCESS_DATA_TYPE_DEPENDENCIES = \
    {'StoichiometricData': ['_stoichiometry', 'lower_bound', 'upper_bound',
                            'subreactions'],

     'ComplexData': ['stoichiometry', 'complex_id', 'modifications'],

     'TranscriptionData': ['modifications', 'subreactions',
                           'nucleotide_sequence', 'RNA_products',
                           'RNA_polymerase', 'rho_dependent'],

     'TranslationData': ['modifications', 'subreactions',
                         'nucleotide_sequence', 'mRNA', 'protein'],

     'tRNAData': ['modifications', 'codon', 'RNA', 'amino_acid',
                  'synthetase', 'synthetase_keff'],

     'TranslocationData': ['enzyme_dict', 'stoichiometry', 'keff',
                           'length_dependent_energy'],

     'PostTranslationData': ['processed_protein_id', 'unprocessed_protein_id',
                             'propensity_scaling', 'aggregation_propensity',
                             'translocation', 'modifications', 'subreactions',
                             'surface_area', 'keq_folding', 'k_folding'],

     'ModificationData': ['stoichiometry', 'enzyme', 'keff'],

     'SubreactionData': ['stoichiometry', 'enzyme', 'keff'],

     'GenericData': ['component_list']
     }

_REQUIRED_METABOLITE_ATTRIBUTES = {"id", "name", "formula"}

_OPTIONAL_METABOLITE_ATTRIBUTES = {"charge", "formula", "compartment",
                                   "_bound", "_constraint_sense"}

# Some metabolite types require additional attributes
_METABOLITE_TYPE_DEPENDENCIES = \
    {'TranscribedGene': ['left_pos', 'right_pos', 'strand', 'RNA_type',
                         'nucleotide_sequence'],
     'ProcessedProtein': ['unprocessed_protein_id']
     }


def get_schema():
    with open(os.path.join(cur_dir, 'JSONSCHEMA'), 'r') as f:
        return json.load(f)


def get_process_data_class(name):
    # Map name in model dictlist to process data class name
    dictlist_to_class = {'stoichiometric_data': 'StoichiometricData',
                         'complex_data': 'ComplexData',
                         'transcription_data': 'TranscriptionData',
                         'translation_data': 'TranslationData',
                         'tRNA_data': 'tRNAData',
                         'translocation_data': 'TranslocationData',
                         'posttranslation_data': 'PostTranslationData',
                         'modification_data': 'ModificationData',
                         'subreaction_data': 'SubreactionData',
                         'generic_data': 'GenericData'}

    try:
        return dictlist_to_class[name]
    except KeyError:
        raise TypeError('Not a valid process_data dictlist')


def _fix_type(value, kind=''):
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
    #if value is None:
    #    return ''
    return value


def _reaction_to_dict(reaction):
    new_reaction = {key: _fix_type(getattr(reaction, key), 'reaction')
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

    special_list = ['_stoichiometry', 'subreactions', 'modifications',
                    'stoichiometry', 'enzyme_dict', 'translocation',
                    'surface_area', 'keq_folding', 'k_folding']

    for attribute in _PROCESS_DATA_TYPE_DEPENDENCIES[process_data_type]:
        if attribute not in special_list:
            data_attribute = getattr(data, attribute)

            new_data[attribute] = _fix_type(data_attribute, 'reaction')

        elif attribute == 'enzyme_dict':
            new_data[attribute] = {}
            for cplx, values in getattr(data, attribute).items():
                new_data[attribute][cplx] = {}
                for property, value in values.items():
                    new_data[attribute][cplx][property] = _fix_type(value)
        else:
            new_data[attribute] = {}
            for metabolite, coefficient in getattr(data, attribute).items():
                new_data[attribute][metabolite] = _fix_type(coefficient)

    return new_data


def _metabolite_to_dict(metabolite):

    # TODO this is incrrectly in model as transcribed gene
    if metabolite.id == 'RNA_degradosome':
        metabolite.left_pos = 1
        metabolite.right_pos = 1
        metabolite.strand = '+'
        metabolite.RNA_type = 'mRNA'

    metabolite_type = metabolite.__class__.__name__

    new_metabolite = {key: _fix_type(getattr(metabolite, key))
                      for key in _REQUIRED_METABOLITE_ATTRIBUTES}

    # Som metabolites require additional information to construct working
    # ME-model
    new_metabolite['metabolite_type'] = {}
    new_metabolite['metabolite_type'][metabolite_type] = {}
    for attribute in _METABOLITE_TYPE_DEPENDENCIES.get(metabolite_type, []):
        metabolite_attribute = getattr(metabolite, attribute)
        new_metabolite['metabolite_type'][metabolite_type][attribute] = \
            metabolite_attribute

    return new_metabolite


def get_attribute_array(dictlist, type):
    if type == 'reaction':
        return [_reaction_to_dict(reaction) for reaction in dictlist]
    elif type == 'process_data':
        return [_process_data_to_dict(data) for data in dictlist]
    elif type == 'metabolite':
        return [_metabolite_to_dict(metabolite) for metabolite in dictlist]
    else:
        raise TypeError('Type must be reaction, process_data or metabolite')


def get_global_info_dict(global_info):
    new_global_info = {}
    for key, value in global_info.items():
        if type(value) != dict:
            new_global_info[key] = _fix_type(value)
        else:
            new_global_info[key] = value
    return new_global_info


def _to_dict(model):

    obj = dict(
        reactions=get_attribute_array(model.reactions, 'reaction'),
        stoichiometric_data=get_attribute_array(model.stoichiometric_data,
                                                'process_data'),
        complex_data=get_attribute_array(model.complex_data, 'process_data'),
        transcription_data=get_attribute_array(model.transcription_data,
                                               'process_data'),
        translation_data=get_attribute_array(model.translation_data,
                                             'process_data'),
        tRNA_data=get_attribute_array(model.tRNA_data, 'process_data'),
        translocation_data=get_attribute_array(model.translocation_data,
                                               'process_data'),
        posttranslation_data=get_attribute_array(model.posttranslation_data,
                                                 'process_data'),
        modification_data=get_attribute_array(model.modification_data,
                                              'process_data'),
        subreaction_data=get_attribute_array(model.subreaction_data,
                                             'process_data'),
        generic_data=get_attribute_array(model.generic_data, 'process_data'),
        metabolites=get_attribute_array(model.metabolites, 'metabolite'),
        global_info=get_global_info_dict(model.global_info)
    )

    return obj


def save_full_me_model_json(model, file_name):
    """
    Save a full JSON version of the ME-model. Saving/loading a model in this
    format can then be loaded to return a ME-model identical to the one saved.

    :param :class:`~cobrame.core.MEModel.MEModel` model:
        A full ME-model

    :param str or file-like object file_name:
        Filename of the JSON output

    :returns JSON-object:
        Full JSON representation of full ME-model
    """

    should_close = False
    if isinstance(file_name, string_types):
        file_name = open(file_name, 'w')
        should_close = True

    json.dump(_to_dict(model), file_name)

    if should_close:
        file_name.close()


def add_metabolite_from_dict(model, metabolite_info):
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
    if metabolite_type != 'ProcessedProtein':
        metabolite_obj = \
            getattr(cobrame, metabolite_type)(metabolite_info['id'])
    else:
        unprocessed_id = \
            metabolite_type_dict['ProcessedProtein']['unprocessed_protein_id']
        metabolite_obj = \
            getattr(cobrame, metabolite_type)(metabolite_info['id'],
                                              unprocessed_id)

    for attribute in _REQUIRED_METABOLITE_ATTRIBUTES:
        setattr(metabolite_obj, attribute, metabolite_info[attribute])

    for attribute in _METABOLITE_TYPE_DEPENDENCIES.get(metabolite_type, []):
        value = metabolite_type_dict[metabolite_type][attribute]
        setattr(metabolite_obj, attribute, value)

    model.add_metabolites([metabolite_obj])


def add_process_data_from_dict(model, process_data_type, process_data_info):
    """
    Builds process_data instances defined in dictionary, then add it to the
    ME-model being constructed.

    Most classes of process_data only require an id and model to initiate them,
    but TranslationData, tRNAData, PostTranslationData and GenericData require
    additional inputs.

    """

    # Create process data instances. Handel certain types individually
    id = process_data_info['id']
    if process_data_type == 'TranslationData':
        mRNA = process_data_info['mRNA']
        protein = process_data_info['protein']
        process_data = \
            getattr(cobrame, process_data_type)(id, model, mRNA, protein)
    elif process_data_type == 'tRNAData':
        amino_acid = process_data_info['amino_acid']
        RNA = process_data_info['RNA']
        codon = process_data_info['codon']
        process_data = \
            getattr(cobrame, process_data_type)(id, model, amino_acid, RNA,
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
        setattr(process_data, attribute, process_data_info[attribute])

    # Some attributes depend on process data type. Set those here.
    for attribute in _PROCESS_DATA_TYPE_DEPENDENCIES.get(process_data_type,
                                                         []):
        value = process_data_info[attribute]
        if type(value) != dict:
            setattr(process_data, attribute, value)
        else:
            for k, v in value.items():
                getattr(process_data, attribute)[k] = v


def add_reaction_from_dict(model, reaction_info):
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
        model.add_reaction(reaction_obj)
    except Exception:
        reaction_obj = model.reactions.get_by_id(reaction_obj.id)
        if reaction_type not in ['SummaryVariable',
                                 'GenericFormationReaction']:
            warn('Reaction (%s) already in model' % reaction_obj.id)

    # These reactions types do not have update functions and need their
    # stoichiometries set explicitly .
    if reaction_type in ['SummaryVariable', 'Reaction']:
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


def full_me_model_from_dict(obj):
    """
    Validate and load JSON representation of the ME-model. This will return
    a full :class:`~cobrame.core.MEModel.MEModel` object identical to the
    one saved.

    :param str or file-like object obj:
        JSON-serialized ME-model

    :returns :class:`~cobrame.core.MEModel.MEModel`:
        Full COBRAme ME-model
    """

    try:
        validate(obj, get_schema())
    except ValidationError:
        raise Exception('Must pass valid ME-model json file')

    model = cobrame.MEModel()

    for k, v in iteritems(obj):
        if k in {'id', 'name', 'global_info'}:
            setattr(model, k, v)

    for metabolite in obj['metabolites']:
        add_metabolite_from_dict(model, metabolite)

    for key in obj:
        if 'data' in key:
            for data in obj[key]:
                add_process_data_from_dict(model, get_process_data_class(key),
                                           data)

    for reaction in obj['reactions']:
        add_reaction_from_dict(model, reaction)

    model.update()

    return model


def load_full_me_model_json(file_name):

    with open(file_name, 'r') as f:
        model_dict = json.load(f)

    return full_me_model_from_dict(model_dict)
