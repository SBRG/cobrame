# COBRAme
COBRA Toolbox for constructing and manipulating genome-scale models of metabolism and expression (ME-models)

## Installation
1. clone the repository
2. run ```python setup.py develop --user```
3. Download the [latest soplex release](https://github.com/SBRG/soplex_cython/releases)
4. Run ```pip install <PATH_TO_SOPLEX_WHL> --user```

## Build ME-model of *E. coli* K-12 MG1655
See [ecoliME](https://github.com/SBRG/ecoliME)

## Requirments
Currently, COBRAme is only supported for Python 2.7, but updating the code for Python 3 is a top priority.

Linux is recommended, with a relatively recent glibc. Ubuntu 14.04 or later
should work. Mac OS X probably works too. Windows has worked in the past, but
is not explicitly supported at this time.

