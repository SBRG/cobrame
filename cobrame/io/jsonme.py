from cobra.io.json import json_schema, save_json_model, load_json_model
from sympy import Basic, sympify
import copy


def save_json_me(me0, file_name, pretty=False):
    """
    Save ME model as json

    model : :class:`~cobrame.core.MEModel.MEmodel` object

    file_name : str or file-like object
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


def load_json_me(file_name):
    """
    Load ME model from json

    file_name : str or file-like object
    """
    me = load_json_model(file_name)

    # Re-convert stoichiometries back to sympy
    for rxn in me.reactions:
        for met in rxn.metabolites:
            s =rxn._metabolites[met]
            try:
                rxn._metabolites[met] = float(s)
            except ValueError:
                rxn._metabolites[met] = sympify(s)
        try:
            rxn.lower_bound = float(rxn.lower_bound)
        except ValueError:
            rxn.lower_bound = sympify(rxn.lower_bound)
        try:
            rxn.upper_bound = float(rxn.upper_bound)
        except ValueError:
            rxn.upper_bound = sympify(rxn.upper_bound)

    for met in me.metabolites:
        b = met._bound
        try:
            met._bound = float(b)
        except ValueError:
            met._bound = sympify(b)

    return me


"""
Add ME objects to the cobrapy json-schema
Particularly, translation_data, etc.
"""
#json_schema
