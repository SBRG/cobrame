from time import time

from cobra.solvers import soplex

from minime.solve.symbolic import *

try:
    from IPython.utils.coloransi import TermColors
    Red = TermColors.Red
    Green = TermColors.Green
    Normal = TermColors.Normal
except ImportError:
    Red = Green = Normal = ""


def batch_growth(me_model, max_mu=2, mu_accuracy=1e-9, verbose=True,
                 **solver_args):
    """Computes maximum growth rate

    This function is appropriate for "batch" growth where nutrient uptake is
    unbounded.

    max_mu: A guess for a growth rate which will be infeasible
    """
    lp = soplex.create_problem(me_model)
    for name, value in iteritems(solver_args):
        lp.set_parameter(name, value)
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
                print "%s%.13f%s\t%s" % (Green, mu, Normal, str(objective))
            if objective > 0:
                feasible_mu.append(mu)
                objectives.append(objective)
                return True
            else:
                assert False
        else:
            infeasible_mu.append(mu)
            if verbose:
                print "%s%.13f%s" % (Red, mu, Normal)
            return False

    start = time()
    # try the edges of binary search
    if not try_mu(0):
        print "0 needs to be feasible"
    while try_mu(max_mu):  # If max_mu was feasible, keep increasing until not
        max_mu += 1
    while infeasible_mu[-1] - feasible_mu[-1] > mu_accuracy:
        try_mu((infeasible_mu[-1] + feasible_mu[-1]) * 0.5)
    if verbose:
        print "completed in %.1f seconds and %d iterations" % \
            (time() - start, len(feasible_mu) + len(infeasible_mu))
    try_mu(feasible_mu[-1])
    me_model.solution = soplex.format_solution(lp, me_model)
    me_model.solution.f = feasible_mu[-1]
