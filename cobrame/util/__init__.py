from __future__ import absolute_import

from sympy import Symbol

# If cobrapy and optlang cannot handle symbolic parameters, assume using
# cobrapy versions <= 0.5.11
try:
    from optlang.interface import SymbolicParameter
except ImportError:
    mu = Symbol("mu", positive=True)
else:
    mu = SymbolicParameter("mu", value=1.)
    mu._assumptions._tell('positive', True)
    mu._assumptions._tell('nonzero', True)
    mu._assumptions._tell('negative', False)
    mu._assumptions._tell('zero', False)
    mu._assumptions._tell('real', True)
    # mu._assumptions._tell('nonnegative', True)
    mu._assumptions['uuid'] = None
