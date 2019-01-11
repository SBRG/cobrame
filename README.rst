|Build_Status| |License| |Documentation| |Docker|

COBRAme
=======

A COBRApy_ extension for constructing and simulating genome-scale models of metabolism and gene expression (ME-models)

For more information on ME-models and the COBRAme ME-model architecture, see the COBRAme readthedocs_.

If using COBRAme or iJL1678b-ME in a publication, please cite: `doi:10.1371/journal.pcbi.1006302 <https://doi.org/10.1371/journal.pcbi.1006302>`_. 

COBRAme with Docker
-------------------
Docker images are available on DockerHub_ which contain pre-installed versions COBRAme, the solvers, all dependencies using Python 3.6. Dockerfiles_ are also available to build Docker containers locally.

The DockerHub images contain a precompiled version the qMINOS solver and SoPlex can additionally be installed if a Docker container is build locally.

Installation
------------

1. clone the repository
2. run ``python setup.py develop --user``

Build ME-model
--------------
A full COBRAme ME-model currently only exists for *E. coli* K-12 MG1655. See the ECOLIme_ package for the data and scripts required to build iJL1678b-ME. After installing ECOLIme_, iJL1678b-ME can be constructed from scratch by running `build_ME_model.py <https://github.com/SBRG/ecolime/tree/master/ecolime>`_.

**Note:** ECOLIme and COBRAme are versioned together (i.e. version 0.0.8 of ECOLIme should be used with version 0.0.8 of COBRAme, etc.)

Solving
-------
ME-models are inherently ill-scaled and thus require extended precision solvers in to be solved accurately. COBRAme currently supports the use of two such solvers: **qMINOS** (128-bit precision) and **SoPlex** (80-bit precision). Both solvers are freely available for academic use. For further examples, refer to `solve_demo.ipynb <https://github.com/SBRG/ecolime/tree/master/ecolime>`_

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

- Python versions 2.7+/3.5+
- COBRApy versions <= 0.5.11. We recommend using 0.5.11.
- Linux is recommended, with a relatively recent glibc. Mac OS X is also supported. Windows has worked in the past, but is not explicitly supported at this time.

.. _readthedocs: http://cobrame.readthedocs.io/
.. _ECOLIme: https://github.com/SBRG/ECOLIme
.. _ZIB: http://soplex.zib.de/
.. _soplex_cython: https://github.com/SBRG/soplex_cython
.. _solvemepy: https://github.com/SBRG/solvemepy
.. _COBRApy: https://github.com/opencobra/cobrapy
.. _DockerFiles: https://github.com/SBRG/cobrame/tree/master/docker
.. _DockerHub: https://hub.docker.com/r/sbrg/cobrame/
.. |Build_Status| image:: https://travis-ci.org/SBRG/cobrame.svg?branch=master
    :target: https://travis-ci.org/SBRG/cobrame
.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://github.com/SBRG/cobrame/blob/master/LICENSE
.. |Documentation| image:: https://readthedocs.org/projects/cobrame/badge/?version=master
    :target: http://cobrame.readthedocs.io/en/master/?badge=master
.. |Docker| image:: https://img.shields.io/docker/build/sbrg/cobrame.svg
    :target: https://hub.docker.com/r/sbrg/cobrame/builds/
