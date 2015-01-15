from sympy import Basic, lambdify
from six import iteritems

from minime import mu


def _compile(expr):
    """compiles a sympy expression, and returns None otherwise"""
    return lambdify(mu, expr) if isinstance(expr, Basic) else None


def compile_expressions(me_model):
    """compiles symbolic exressions of mu to functions

    The compiled expressions dict has the following key value pairs:
    (met_index, rxn_index): stoichiometry,
    (None, rxn_index): (lower_bound or None, upper_bound or None)
    (met_index, None): met_bound

    Any non-None value in the dict can be assumed to be a compiled function
    """
    expressions = {}
    for i, r in enumerate(me_model.reactions):
        # stoichiometry
        for met, stoic in iteritems(r._metabolites):
            if isinstance(stoic, Basic):
                expressions[(me_model.metabolites.index(met), i)] = \
                    lambdify(mu, stoic)
        bounds = (_compile(r.lower_bound), _compile(r.upper_bound))
        if bounds[0] is not None or bounds[1] is not None:
            expressions[(None, i)] = bounds
    for i, metabolite in enumerate(me_model.metabolites):
        compiled = _compile(metabolite._bound)
        if compiled is not None:
            expressions[(i, None)] = compiled
    return expressions


def substitute_mu(lp, mu, compiled_exressions):
    """substitute mu into a constructed LP"""
    for index, expr in iteritems(compiled_exressions):
        if index[0] is None:  # reaction bounds
            print "not implemented"
        elif index[1] is None:  # metabolite _bound
            print "not implemented"
        else:  # stoichiometry
            lp.change_coefficient(index[0], index[1], expr(mu))
