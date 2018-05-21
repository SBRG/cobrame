COBRAme, ECOLIme, and qMINOS using Docker
=========================================

Builds of the latest versions of COBRAme and ECOLIme_ can be found on Docker Hub. This build also includes a compiled version of the qMINOS solver that can be used to solve ME-models using solvemepy_.

This container will initialize  which contains a json and pickle of

Docker is cross-platform and allows Windows, Mac and Linux users to run ME-models without going through the process of installing each solver and dependency locally.


Installation
------------
To run a Docker container with everything required to build and solve the iJL1678b model of *E. coli* K-12 MG1655 using qMINOS:

1. Install Docker (https://docs.docker.com/install/)
2. In the command line, run ``docker run -p 8888:8888 --rm -i -v $(pwd):/mount_point/ -t sbrg/cobrame:master bash``.

This will initiate a Docker container (virtual machine) into the ``/home/meuser`` directory and mount the contents of the directory where the command was ran into the docker container at ``/mount_point/``.

3. To start a jupyter notebook, run ``jupyter notebook --ip=0.0.0.0 --port=8888``. Point the browser to ``localhost:8888`` and input the provided token to access notebook.


Building Docker locally w/ SoPlex
---------------------------------
To build a Docker image that can solve ME-models with SoPlex as well as qMINOS

1. Install Docker (https://docs.docker.com/install/)
2. Download SoPlex version 3.1.1 (http://soplex.zib.de/#download)
3. ``cp soplex-3.1.1.tgz [cobrame root]/docker``
4. ``cd [cobrame root]/docker``
5. ``docker build -t [repository:tag name] .``

where the ``[repository:tag name]`` can be decided by the user

6. To run the built Docker image, ``docker run -p 8888:8888 --rm -i -v $(pwd):/mount_point/ -t [repository:tag name] bash``

This will initialize a Docker container with everything required to use ME-models and to solve them using qMINOS and SoPlex.

.. _ECOLIme: https://github.com/SBRG/ECOLIme
.. _ZIB: http://soplex.zib.de/
.. _soplex_cython: https://github.com/SBRG/soplex_cython
.. _solvemepy: https://github.com/SBRG/solvemepy
.. _COBRApy: https://github.com/opencobra/cobrapy