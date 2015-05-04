from sympy import Basic, lambdify
from six import iteritems

from minime import mu


def _compile(expr):
    """compiles a sympy expression"""
    return lambdify(mu, expr) if isinstance(expr, Basic) else expr


def _eval(expr, mu):
    """evaluate the expression as a function of mu

    will return the expr if it is not a function"""
    return expr(mu) if callable(expr) else expr


def compile_expressions(me_model):
    """compiles symbolic exressions of mu to functions

    The compiled expressions dict has the following key value pairs:
    (met_index, rxn_index): stoichiometry,
    (None, rxn_index): (lower_bound, upper_bound)
    (met_index, None): (met_bound, met_constraint_sense)

    """
    expressions = {}
    for i, r in enumerate(me_model.reactions):
        # stoichiometry
        for met, stoic in iteritems(r._metabolites):
            if isinstance(stoic, Basic):
                expressions[(me_model.metabolites.index(met), i)] = \
                    lambdify(mu, stoic)
        # If either the lower or upper reaction bounds are symbolic
        if isinstance(r.lower_bound, Basic) or \
                isinstance(r.upper_bound, Basic):
            expressions[(None, i)] = (_compile(r.lower_bound),
                                      _compile(r.upper_bound))
    # Metabolite bound
    for i, metabolite in enumerate(me_model.metabolites):
        if isinstance(metabolite._bound, Basic):
            expressions[(i, None)] = (_compile(metabolite._bound),
                                      metabolite._constraint_sense)
    return expressions


def substitute_mu(lp, mu, compiled_exressions):
    """substitute mu into a constructed LP"""
    for index, expr in iteritems(compiled_exressions):
        if index[0] is None:  # reaction bounds
            lp.change_variable_bounds(index[1],
                                      _eval(expr[0], mu), _eval(expr[1], mu))
        elif index[1] is None:  # metabolite _bound
            lp.change_constraint(index[0], expr[1], _eval(expr[0], mu))
        else:  # stoichiometry
            lp.change_coefficient(index[0], index[1], _eval(expr, mu))
