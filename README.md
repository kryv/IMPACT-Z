IMPACT: Integrated Map and Particle Accelerator Tracking Code
=============================================================

This branch contains the original IMPACT source code projected by
Ji Qiang from LBNL with modifications to support development of FRIB.

## Pre-requisites ##
Needs mpi library for fortran source, also scipy and mpi4py for Python interface.
The nosetests test runner is used for 'make test' if present.

For Python2.x,

```sh
$ apt-get install libopenmpi-dev python-scipy python-nose python-mpi4py cmake
```

For Python3.x,

```sh
$ apt-get install libopenmpi-dev python3-scipy python3-nose python3-mpi4py cmake
```

## Installation ##
`git clone` this repository, and

```sh
$ cd impact
$ mkdir build
$ cd build
$ cmake ..
$ make
```

### Main contents of `build` directory ###
  * bin/Impact - Original IMPACT executable file
  * bin/AdvImpact - Advanced IMPACT executable file
  * lib/ - Library files for IMPACT
  * modules/ - Fortran modules directory
  * python/ - python module for IMPACT

## Documentation ##
* Advanced IMPACT usage ([html](Docs/html/), [pdf](Docs/AdvIMPACT.pdf))
* [Original documentation](Docs/ImpactZv1Readme3.txt)
