# External libraries

The libraries described in the following provide additional functionality
but are not necessary to install and work DuMux. If you are going to use an external library,
check the information provided on the [DUNE website](https://www.dune-project.org/doc/external-libraries/).
The DUNE website provides information for [modules providing grid implementations](https://www.dune-project.org/groups/grid/)
or other [extension modules](https://www.dune-project.org/groups/extension/).

Installing an external library can require additional libraries which are also used by DUNE.
For some libraries, such as BLAS or MPI, multiple versions can be installed on the system.
Make sure that it uses the same library as DUNE when configuring the external library.

Some of the libraries are then compiled within that directory and are not installed in a different place,
but DUNE may need to know their location. See [installation instructions](#install-external-dependencies-via-script)
for tips on how to install external libraries.

In the following grouped lists, you can find some external modules and external libraries.

[TOC]

## Grid managers

* `dune-alugrid`: Grid library, comes as a DUNE module. The parallel version needs also a graph partitioner, such as ParMETIS.
Download: https://gitlab.dune-project.org/extensions/dune-alugrid.git
* `dune-foamgrid`: External grid module. One- and two-dimensional grids in a physical space of arbitrary dimension; non-manifold grids, growth, element paramterizations, and movable vertices. This makes FoamGrid the grid data structure of choice for simulating structures such as foams, discrete fracture networks, or network flow problems. Download: https://gitlab.dune-project.org/extensions/dune-foamgrid.git
* `opm-grid`: opm-grid is a DUNE module supporting grids in a corner-point format. Download: https://github.com/OPM/opm-grid.git and https://github.com/OPM/opm-common.git
* `dune-subgrid`: The dune-subgrid module is a meta-grid implementation that allows to mark elements of another hierarchical dune grid and use this sub-grid just like a regular grid. The set of marked elements can then be accessed as a hierarchical dune grid in its own right. Dune-Subgrid provides the full grid interface including adaptive mesh refinement. Download: https://git.imp.fu-berlin.de/agnumpde/dune-subgrid.git
* `dune-spgrid`: The DUNE module dune-spgrid provides a structured, parallel grid and supports periodic boundary conditions.
Download: https://gitlab.dune-project.org/extensions/dune-spgrid.git
* `dune-uggrid`: External library for use as grid. UG is a toolbox for unstructured grids, released under GPL. To build UG the tools lex/yacc or the GNU variants of flex/bison must be provided. Download: https://gitlab.dune-project.org/staging/dune-uggrid.git

## Solver libraries

* `SuperLU`: External library for solving linear equations. SuperLU is a general purpose library for the direct solution of large, sparse, non-symmetric systems of linear equations. Download: http://crd.lbl.gov/~xiaoye/SuperLU
* `UMFPack`: External library for solving linear equations. It is part of `SuiteSparse`. See: http://faculty.cse.tamu.edu/davis/suitesparse.html. On Debian/Ubuntu you can install the package `libsuitesparse-dev`.

## Parallel computing (distributed memory)

* `MPI`: The distributed parallel version of DUNE and also some of the external dependencies need MPI when they are going to be built for parallel computing. `OpenMPI` and `MPICH` in a recent version have been reported to work.

## Parallel computing (shared memory)

If one of the following libraries is installation, multi-threaded code is enabled in DuMux:

* `OpenMP`: OpenMP usually comes with the compiler
* `OneTBB`: TBB, see https://github.com/oneapi-src/oneTBB
* `C++ parallel algorithms`: see https://en.cppreference.com/w/cpp/algorithm
* `Kokkos`: see https://github.com/kokkos/kokkos

The chosen backend can be set by passing the variable `DUMUX_MULTITHREADING_BACKEND` to CMake,
or by modifying it in the CMake cache (in `dumux/build-cmake` run `ccmake .`).
When running `dunecontrol` or CMake, the output contains a line like

```
-- Dumux multithreading backed: TBB
```

which in this case shows that the selected backend is TBB.

## Linear algebra libraries

The following are dependencies of some of the used libraries. You will need them depending on which modules of DUNE and which external libraries you use.
* `BLAS`: `SuperLU` makes use of `BLAS`. Thus install OpenBLAS, GotoBLAS2, ATLAS, non-optimized BLAS or BLAS provided by a chip manufacturer.
Take care that the installation scripts select the intended version of BLAS.
* `METIS` and `ParMETIS`: This are dependencies of ALUGrid and can be used with UG, if run in parallel.
