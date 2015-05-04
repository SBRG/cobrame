from time import time
from warnings import warn

from cobra.solvers import solver_dict
from cobra.flux_analysis.variability import calculate_lp_variability
from six import string_types, iteritems

try:
    import soplex
except ImportError as e:
    soplex = None
    warn("soplex import failed with error '%s'" % e.message)

from minime.solve.symbolic import *

try:
    from IPython.utils.coloransi import TermColors
    Red = TermColors.Red
    Green = TermColors.Green
    Normal = TermColors.Normal
except ImportError:
    Red = Green = Normal = ""


def get_ME_solver(solver=None):
    if solver is None:
        if soplex is None:
            raise RuntimeError("soplex not installed")
        return soplex
    elif isinstance(solver, string_types):
        return solver_dict[solver]
    else:
        return solver


def binary_search(me_model, min_mu=0, max_mu=2, mu_accuracy=1e-9,
                  solver=None, verbose=True, compiled_expressions=None,
                  **solver_args):
    """Computes maximum feasible growth rate (mu) through a binary search

    The objective function of the model should be set to a dummy
    reaction which forces translation of a dummy protein.

    max_mu: A guess for a growth rate which will be infeasible
    min_mu: A guess for a growth rate which will be feasible
    mu_accuracy: The final error in mu after the binary search
    verbose: will print out each mu in the binary search
    compiled_expressions: precompiled symbolic expressions in the model

    """
    solver = get_ME_solver(solver)
    for name, value in iteritems(solver_args):
        lp.set_parameter(name, value)
    if compiled_expressions is None:
        compiled_expressions = compile_expressions(me_model)
    feasible_mu = []
    infeasible_mu = []
    objectives = []

    def try_mu(mu):
        substitute_mu(lp, mu, compiled_expressions)
        lp.solve_problem()
        status = lp.get_status()
        if status == "optimal":
            objective = lp.get_objective_value()
            if verbose:
                print("%s%.13f%s\t%s" % (Green, mu, Normal, str(objective)))
            if objective > 0:
                feasible_mu.append(mu)
                objectives.append(objective)
                return True
            else:
                assert False
        else:
            infeasible_mu.append(mu)
            if verbose:
                print("%s%.13f%s" % (Red, mu, Normal))
            return False

    start = time()
    # try the edges of binary search
    if not try_mu(min_mu):
        # Try 0 if min_mu failed
        if min_mu == 0 or not try_mu(0):
            raise ValueError("0 needs to be feasible")
    while try_mu(max_mu):  # If max_mu was feasible, keep increasing until not
        max_mu += 1
    while infeasible_mu[-1] - feasible_mu[-1] > mu_accuracy:
        try_mu((infeasible_mu[-1] + feasible_mu[-1]) * 0.5)
    try_mu(feasible_mu[-1])
    me_model.solution = solver.format_solution(lp, me_model)
    me_model.solution.f = feasible_mu[-1]
    if verbose:
        print("completed in %.1f seconds and %d iterations" %
              (time() - start, len(feasible_mu) + len(infeasible_mu)))
    return me_model.solution


def create_lP_at_growth_rate(me_model, growth_rate, compiled_expressions=None,
                             solver=None, **solver_args):
    solver = get_ME_solver(solver)
    lp = solver.create_problem(me_model)
    for name, value in iteritems(solver_args):
        lp.set_parameter(name, value)
    # substitute in values
    if compiled_expressions is None:
        compiled_expressions = compile_expressions(me_model)
    substitute_mu(lp, growth_rate, compiled_expressions)
    return (lp, solver)


def solve_at_growth_rate(me_model, growth_rate, **solver_args):
    lp, solver = create_lP_at_growth_rate(me_model, growth_rate,
                                          **solver_args)
    # solve and return
    solver.solve_problem(lp)
    me_model.solution = solver.format_solution(lp, me_model)
    if me_model.solution.status == "optimal":
        me_model.solution.f = growth_rate
    return me_model.solution


def fva(me_model, growth_rate, reaction_list, **solver_args):
    lp, solver = create_lP_at_growth_rate(me_model, growth_rate,
                                          **solver_args)
    return calculate_lp_variability(lp, solver, me_model, reaction_list)
