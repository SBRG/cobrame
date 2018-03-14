COBRAme, ECOLIme, and qMINOS using Docker
=========================================

Builds of the latest versions of COBRAme and ECOLIme_ can be found on Docker Hub. This build also includes a compiled version of the qMINOS solver that can be used to solve ME-models using solvemepy_.

Docker is cross-platform and allows Windows, Mac and Linux users to run ME-models without going through the process of installing each solver and dependency locally.


Installation
------------
To run a Docker container with everything required to build and solve the iJL1678b model of *E. coli* K-12 MG1655 using qMINOS:

1. Install Docker (https://docs.docker.com/install/)
2. In the command line, run ``docker run --rm -i -v $(pwd):/home/meuser/ -t sbrg/cobrame:everything bash``.

This will initiate a Docker container (virtual machine) and copy everything in the directory where the command was run into the docker container at ``/workdir/``

3. To start a jupyter notebook, run ``sh /home/meuser/run_jupyter.sh``


Building Docker locally w/ SoPlex
---------------------------------
To build a Docker image that can solve ME-models with SoPlex as well as qMINOS

1. Install Docker (https://docs.docker.com/install/)
2. Download SoPlex version 3.1.1 (http://soplex.zib.de/#download)
3. ``cp soplex-3.1.1.tgz [cobrame root]/docker``
4. ``cd [cobrame root]/docker``
5. ``docker build -t [repository:tag name] .``

where the ``[repository:tag name]`` can be decided by the user

6. To run the built Docker image, ``docker run --rm -i -v $(pwd):/home/meuser/ -t [repository:tag name] bash``

This will initialize a Docker container with everything required to use ME-models and to solve them using qMINOS and SoPlex.

.. _ECOLIme: https://github.com/SBRG/ECOLIme
.. _ZIB: http://soplex.zib.de/
.. _soplex_cython: https://github.com/SBRG/soplex_cython
.. _solvemepy: https://github.com/SBRG/solvemepy
.. _COBRApy: https://github.com/opencobra/cobrapy