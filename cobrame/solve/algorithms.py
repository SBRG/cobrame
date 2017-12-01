from __future__ import print_function, division, absolute_import

from math import log
from os.path import join
from tempfile import mkdtemp
from time import time
from warnings import warn

from cobra.flux_analysis.variability import calculate_lp_variability
from cobra.solvers import solver_dict
from six import string_types, iteritems

try:
    import soplex
except ImportError as e:
    soplex = None
    warn("soplex import failed with error '%s'" % e)

from cobrame.solve.symbolic import *

try:
    from IPython.utils.coloransi import TermColors
    Red = TermColors.Red
    Green = TermColors.Green
    Normal = TermColors.Normal
except ImportError:
    Red = Green = Normal = ""


def get_me_solver(solver=None):
    if solver is None or solver == 'soplex':
        if soplex is None:
            raise RuntimeError("soplex not installed")
        return soplex
    elif isinstance(solver, string_types):
        return solver_dict[solver]
    else:
        return solver


def binary_search(me_model, min_mu=0, max_mu=2, mu_accuracy=1e-9,
                  solver=None, verbose=True, compiled_expressions=None,
                  debug=True, reset_obj=False, **solver_args):
    """Computes maximum feasible growth rate (mu) through a binary search

    The objective function of the model should be set to a dummy
    reaction which forces translation of a dummy protein.

    :param float max_mu: A guess for a growth rate which will be infeasible
    :param float min_mu: A guess for a growth rate which will be feasible
    :param float mu_accuracy: The final error in mu after the binary search
    :param boolean verbose: will print out each mu in the binary search
    :param dict compiled_expressions: precompiled symbolic expressions in the
        model

    """

    if solver is not None:
        debug = False  # other solvers can't handle debug mode
    solver = get_me_solver(solver)
    lp = solver.create_problem(me_model)
    # reset the objective for faster feasibility solving
    if reset_obj:
        for i, reaction in enumerate(me_model.reactions):
            if reaction.objective_coefficient != 0:
                solver.change_variable_objective(lp, i, 0)
    for name, value in iteritems(solver_args):
        solver.set_parameter(lp, name, value)
    if compiled_expressions is None:
        compiled_expressions = compile_expressions(me_model)
    feasible_mu = []
    infeasible_mu = []

    # String formatting for display
    str_places = int(abs(round(log(mu_accuracy)/log(10)))) + 1
    num_format = "%." + str(str_places) + "f"
    mu_str = "mu".ljust(str_places + 2)
    if debug:
        verbose = True
        save_dir = mkdtemp()
        print("LP files will be saved in " + save_dir)
        filename_base = join(save_dir, "me_mu_" + num_format + ".lp").encode()
    if verbose:
        success_str_base = Green + num_format + "\t+" + Normal
        failure_str_base = Red + num_format + "\t-" + Normal
        if debug:
            print("%s\tstatus\treset\ttime\titer\tobj" % mu_str)
        else:
            print("%s\tstatus" % mu_str)

    def try_mu(mu):
        if mu == 0 and me_model.global_info.get('k_deg', 0) != 0:
            warn('Due to mRNA degradation constraint formulation the model is '
                 'infeasible at mu = 0. Using mu = .1 instead.')
            mu = .1
        substitute_mu(lp, mu, compiled_expressions, solver)
        if debug:
            lp.write(filename_base % mu, rational=True)
        solver.solve_problem(lp)
        status = solver.get_status(lp)
        if debug:
            reset_basis = lp.reset_basis
            obj = str(lp.get_objective_value()) \
                if status == "optimal" else ""
            debug_str = "\t%s\t%.2f\t%d\t%s" % \
                (reset_basis, lp.solveTime, lp.numIterations, obj)
        else:
            debug_str = ""
        if status == "optimal":
            if verbose:
                print(success_str_base % mu + debug_str)
            feasible_mu.append(mu)
            return True
        else:
            infeasible_mu.append(mu)
            if verbose:
                if status != "infeasible" and debug:
                    debug_str += "(status: %s)" % status
                print(failure_str_base % mu + debug_str)
            return False

    start = time()
    # find highest possible value of mu try the edges of binary search
    if not try_mu(min_mu):
        # Try 0 if min_mu failed
        if min_mu == 0 or not try_mu(0):
            raise ValueError("0 needs to be feasible")
    while try_mu(max_mu):  # If max_mu was feasible, keep increasing
        max_mu += 1
    while infeasible_mu[-1] - feasible_mu[-1] > mu_accuracy:
        try_mu((infeasible_mu[-1] + feasible_mu[-1]) * 0.5)
    # now we want to solve with the objective
    if reset_obj:
        for i, reaction in enumerate(me_model.reactions):
            if reaction.objective_coefficient != 0:
                solver.change_variable_objective(
                    lp, i, reaction.objective_coefficient)
    try_mu(feasible_mu[-1])
    me_model.solution = solver.format_solution(lp, me_model)
    me_model.solution.f = feasible_mu[-1]

    if verbose:
        print("completed in %.1f seconds and %d iterations" %
              (time() - start, len(feasible_mu) + len(infeasible_mu)))
    return me_model.solution


def create_lp_at_growth_rate(me_model, growth_rate, compiled_expressions=None,
                             solver=None, **solver_args):

    if growth_rate == 0 and me_model.global_info.get('k_deg', 0) != 0:
        warn('Due to mRNA degradation constraint formulation the model is '
             'infeasible at mu = 0. Using mu = .1 instead.')
        growth_rate = .1
    solver = get_me_solver(solver)
    lp = solver.create_problem(me_model)
    for name, value in iteritems(solver_args):
        lp.set_parameter(name, value)
    # substitute in values
    if compiled_expressions is None:
        compiled_expressions = compile_expressions(me_model)
    substitute_mu(lp, growth_rate, compiled_expressions, solver)
    return (lp, solver)


def solve_at_growth_rate(me_model, growth_rate, **solver_args):
    lp, solver = create_lp_at_growth_rate(me_model, growth_rate,
                                          **solver_args)
    # solve and return
    solver.solve_problem(lp)
    me_model.solution = solver.format_solution(lp, me_model)
    if me_model.solution.status == "optimal":
        me_model.solution.f = growth_rate
    return me_model.solution


def fva(me_model, growth_rate, reaction_list, skip_check=False, **solver_args):
    # store objective
    if skip_check:
        obj = {}
    else:
        obj = {}
        for r in me_model.reactions:
            if r.objective_coefficient != 0:
                obj[r] = r.objective_coefficient
                r.objective_coefficient = 0

    lp, solver = create_lp_at_growth_rate(me_model, growth_rate,
                                          **solver_args)
    result = calculate_lp_variability(lp, solver, me_model, reaction_list)

    # restore the objective value
    for r, v in iteritems(obj):
        r.objective_coefficient = v

    return result
