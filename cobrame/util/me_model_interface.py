from __future__ import absolute_import, print_function

import six


"""
Optlang Interface for ME-model. Overrides some of the functions in the
high-level interface. If version of cobrapy < 0.6.0, this import does nothing.
"""

try:
    from optlang.interface import SymbolicParameter
except ImportError:
    pass
else:
    from optlang import interface
    from optlang.util import inheritdocstring

    @six.add_metaclass(inheritdocstring)
    class Variable(interface.Variable):
        def __init__(self, name, lb=None, ub=None, type="continuous", *args,
                     **kwargs):
            if type != "continuous":
                raise ValueError("ME-models require continuous variables.")
            super(Variable, self).__init__(name, lb, ub, type, *args, **kwargs)

    @six.add_metaclass(inheritdocstring)
    class Constraint(interface.Constraint):

        def __init__(self, expression, sloppy=False, *args, **kwargs):
            super(Constraint, self).__init__(expression, sloppy=sloppy, *args,
                                             **kwargs)

        def set_linear_coefficients(self, coefficients):
            return

    @six.add_metaclass(inheritdocstring)
    class Objective(interface.Objective):
        def __init__(self, expression, sloppy=False, **kwargs):
            super(Objective, self).__init__(expression, sloppy=sloppy,
                                            **kwargs)

        def set_linear_coefficients(self, coefficients):
            return

    @six.add_metaclass(inheritdocstring)
    class OptimizationExpression(interface.OptimizationExpression):
        def __init__(self, expression, sloppy=False, **kwargs):
            super(OptimizationExpression, self).__init__(expression,
                                                         sloppy=sloppy,
                                                         **kwargs)

        def set_linear_coefficients(self, coefficients):
            return

    @six.add_metaclass(inheritdocstring)
    class Configuration(interface.MathematicalProgrammingConfiguration):
        def __init__(self, *args, **kwargs):
            super(Configuration, self).__init__(*args, **kwargs)

    @six.add_metaclass(inheritdocstring)
    class Model(interface.Model):
        def __init__(self, problem=None, *args, **kwargs):

            super(Model, self).__init__(*args, **kwargs)
