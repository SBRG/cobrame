|Build_Status| |License| |Documentation|

COBRAme
=======

COBRA Toolbox for constructing and manipulating genome-scale models of metabolism and expression (ME-models)

For more information on ME-models and the COBRAme ME-model architecture, see the COBRAme readthedocs_.

If using COBRAme or *i*JL1678b-ME in a publication, please cite: `doi:10.1101/106559 <https://doi.org/10.1101/106559>`_. Note that the model may be edited and updated until it is peer-reviewed and published.

Installation
------------

1. clone the repository
2. run ``python setup.py develop --user``

Build ME-model
--------------
A full COBRAme ME-model currently only exists for *E. coli* K-12 MG1655. See the ECOLIme_ package for the data and scripts required to build *i*JL1678b-ME. After installing ECOLIme_, *i*JL1678b-ME can be constructed from scratch by running `build_ME_model.py <https://github.com/SBRG/ecolime/tree/master/ecolime>`_.

Solving
-------
ME-models are inherently ill-scaled and thus require extended precision solvers in to be solved accurately. COBRAme currently supports the use of two such solvers: **qMINOS** (128-bit precision) and **SoPlex** (80-bit precision). Both solvers are freely avaialbe for academic use. For further examples, refer to `solve_demo.ipynb <https://github.com/SBRG/ecolime/tree/master/ecolime>`_

qMINOS
~~~~~~

To install qMINOS for use with COBRAme:

1. Obtain the qMINOS source code from Prof. Michael A. Saunders at Stanford University
2. Download and install the solvemepy_ extension
3. Once installed, ME-models can be solved using a bisection routine by running the following code:

::

  from qminospy.me1 import ME_NLP1

  me_nlp = ME_NLP1(me_model, growth_key='mu')
  muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-6, mumax=1.5)
  me_model.solution.f = me_model.solution.x_dict['biomass_dilution']
  


SoPlex
~~~~~~

To install SoPlex for use with COBRAme:

1. Download the SoPlex source code from ZIB_
2. Download and install the soplex_cython_ extension 
3. Once installed, ME-models can be solved using a binary search routine by running the following code:

::

  from cobrame.solve.algorithms import binary_search
  
  binary_search(me_model, min_mu=.1, max_mu=1.5, debug=True, mu_accuracy=1e-6)


Requirements
------------

COBRAme and its extensions require:

- Python versions >= 2.7
- Linux is recommended, with a relatively recent glibc. Mac OS X is also supported. Windows has worked in the past, but is not explicitly supported at this time.

.. _readthedocs: http://cobrame.readthedocs.io/en/stable/
.. _ECOLIme: https://github.com/SBRG/ECOLIme
.. _ZIB: http://soplex.zib.de/
.. _soplex_cython: https://github.com/SBRG/soplex_cython
.. _solvemepy: https://github.com/SBRG/solvemepy
.. |Build_Status| image:: https://travis-ci.org/SBRG/cobrame.svg?branch=master
    :target: https://travis-ci.org/SBRG/cobrame
.. |License| image:: https://img.shields.io/github/license/mashape/apistatus.svg
    :target: https://github.com/SBRG/cobrame/blob/master/LICENSE
.. |Documentation| image:: https://readthedocs.org/projects/cobrame/badge/?version=master
    :target: http://cobrame.readthedocs.io/en/master/?badge=master
