msParSm
=======

The _msParSm_ application is an evolution of _msPar_ <sup>[(1)][1]</sup>, the parallel version of the coalescent simulation program 
_ms_ <sup>[(2)][2]</sup>, 
which removes the limitation for simulating long stretches of DNA sequences with large recombination rates, 
without compromising the accuracy of the standard coalescence.

## Pre-requisites
- Linux GNU Compiler 4.9.1 (or greater)
- OpenMPI 1.10.1 (other releases in the branch 1.10 should be fine)
    - Version 1.8.x could potentially be fine, but please notice that _msParSm_ was not fully tested with such version.
- CMake 3.5.1 (or greather) **OR** GNU Make 3.81 (or greater)

## How to Build
There are two ways for building _msParSm_: CMake and Make. If you have installed CMAKE with version greater than 3.5.0,
then go with CMake, otherwise you should use Make.

### CMake
```bash
cmake <src-path> -DCMAKE_INSTALL_PREFIX=<install-path>
```

### Make
```bash
make install
```
Binary files will be put into the `bin` folder (which is already _git ignored_).

## How to Use
Usage is the mostly the same as with traditional _ms_, but you need to run it through _OpenMPI_. Next example
will run the application using 4 threads:

```bash
mpirun -n 4 bin/msparsm 10 20 -seeds 40328 19150 54118 -t 100 -r 100 100000 -I 2 2 8 -eN 0.4 10.01 -eN 1 0.01 -en 0.25 2 0.2 -ej 3 2 1 -T > results.out
```

[1]: http://link.springer.com/chapter/10.1007/978-3-642-54420-0_32
[2]: http://home.uchicago.edu/~rhudson1/popgen356/OxfordSurveysEvolBiol7_1-44.pdf