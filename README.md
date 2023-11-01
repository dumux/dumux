<img src="https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/b9550fddee9b6d6dddb8fcae8e0d14440f7a55eb/doc/logo/dumux_logo_hires_whitebg.png" alt="dumux logo" width="400"/>

# What is DuMux?

[DuMu<sup>x</sup>][0] is a simulation framework with a focus on
finite volume discretization methods, model coupling for multi-physics applications,
and flow and transport applications in porous media.

DuMu<sup>x</sup> is based on the [DUNE][1] framework from which it uses
the versatile grid interface, vector and matrix types, geometry and local basis functions, and linear solvers.
DuMu<sup>x</sup> then provides

* Finite volume discretizations (Tpfa, Mpfa, Staggered) and control-volume finite element discretization schemes
* A flexible system matrix assembler and approximation of the Jacobian matrix by numeric differentation
* A customizable Newton method implementation including line search and various stopping criteria
* Many preimplemented models (Darcy-scale porous media flow, Navier-Stokes, Geomechanics, Pore network models, Shallow water equations) and constitutive models
* A multi-domain framework for model coupling suited to couple subproblems with different discretizations/domains/physics/dimensions/... and create monolithic solvers

DuMu<sup>x</sup> has been applied to model complex and non-linear phenomena,
such as $\mathrm{CO}_2$ sequestration, soil remediation, reactive transport and precipitaion phenomena,
drug delivery in cancer therapy, flow in micro-fluidics, root-soil interaction,
flow in fractured porous media, atmosphere-soil flow interaction, evaporation and more.
Please have a look at our journal publications
(see below: [How to cite](#how-to-cite))
for a more detailed description of the goals the development history
and motivations behind DuMu<sup>x</sup>.

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

[TOC]

# Overview

The following resources are useful to get started with DuMu<sup>x</sup>:

* [Installation guide][3]
* [Getting started guide](https://dumux.org/docs/doxygen/master/getting-started.html)
* [Documentation](https://dumux.org/docs/doxygen/master/),
* [DuMu<sup>x</sup> course materials](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/tree/master),
* [Examples](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples), with detailed description of code and results,
* [Mailing list](https://listserv.uni-stuttgart.de/mailman/listinfo/dumux),
* [Changelog](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CHANGELOG.md), where all changes between different release versions are listed and explained.

Some helpful code snippets are available in the [Wiki](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/wikis/home).

Automated testing of installation: [![installation testing pipeline](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-test-installation/badges/main/pipeline.svg)](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-test-installation/-/pipelines?page=1&scope=all&ref=main)

# License

[![REUSE status](https://api.reuse.software/badge/git.iws.uni-stuttgart.de/dumux-repositories/dumux)](https://api.reuse.software/info/git.iws.uni-stuttgart.de/dumux-repositories/dumux)

DuMu<sup>x</sup> is licensed under the terms and conditions of the GNU General
Public License (GPL) version 3 or - at your option - any later
version. The GPL can be [read online][4] or in the [LICENSE.md](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/afb7f2296d84fd2367c612a1084d9b47ff85a260/LICENSE.md) file
provided in the topmost directory of the DuMu<sup>x</sup> source code tree.

Please note that DuMu<sup>x</sup>' license, unlike DUNE's, does *not* feature a
template exception to the GNU General Public License. This means that
you must publish any source code which uses any of the DuMu<sup>x</sup> header
files if you want to redistribute your program to third parties. If
this is unacceptable to you, please [contact us][5] for a commercial
license.

See the file [LICENSE.md](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/afb7f2296d84fd2367c612a1084d9b47ff85a260/LICENSE.md) for copying permissions.
For a curated list of contributors see [AUTHORS.md](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/afb7f2296d84fd2367c612a1084d9b47ff85a260/AUTHORS.md).
If you notice that a contributor is missing on the list,
please [contact us][5] or open a merge request adding the name.

# How to cite

DuMu<sup>x</sup> is research software and developed at research institutions.
You can cite **specific releases** via [**DaRUS**](https://darus.uni-stuttgart.de/dataverse/iws_lh2_dumux) (from 3.6) or **Zenodo**:
[![zenodo badge](https://zenodo.org/badge/DOI/10.5281/zenodo.2479594.svg)](https://doi.org/10.5281/zenodo.2479594). You can also cite individual code files or even lines via **Software Heritage**: [![sw](https://archive.softwareheritage.org/badge/origin/https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git/)](https://archive.softwareheritage.org/swh:1:dir:e947c9ac369afd90195080e4a06bbde2e1e150ca;origin=https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git;visit=swh:1:snp:3cf49b55de0218903103d428c378e356d7d4d082;anchor=swh:1:rev:11871e4abf619d4cb3f938aedd7a2dea47ce1e87)


If you are using DuMu<sup>x</sup> in scientific publications and in
the academic context, please cite (at least one of)
our publications:

**DuMux 3 – an open-source simulator for solving flow and transport problems in porous media with a focus on model coupling.**
*Computers & Mathematics with Applications*, 81, 423-443, (2021).
[![dumuxCAMWAbadge](https://img.shields.io/badge/DOI-10.1016%2Fj.camwa.2020.02.012-blue)](https://doi.org/10.1016/j.camwa.2020.02.012) [PDF][6]

```bib
@article{Koch2021,
    doi = {10.1016/j.camwa.2020.02.012},
    year = {2021}, volume = {81}, pages = {423--443},
    publisher = {Elsevier {BV}},
    author = {Timo Koch and Dennis Gläser and Kilian Weishaupt and others},
    title = {{DuMux} 3 {\textendash} an open-source simulator for solving flow and transport problems in porous media with a focus on model coupling},
    journal = {Computers \& Mathematics with Applications}}
```

**DuMux: DUNE for multi-{phase,component,scale,physics,…} flow and transport in porous media.**
*Advances in Water Resources*, 34(9), 1102–1112, (2011)
[![dumuxAWRbadge](https://img.shields.io/badge/DOI-10.1016%2Fj.advwatres.2011.03.007-blue)](https://doi.org/10.1016/j.advwatres.2011.03.007) [PDF][2]

```bib
@article{Flemisch2011,
    doi = {10.1016/j.advwatres.2011.03.007},
    year = {2011}, volume = {34}, number = {9}, pages = {1102--1112},
    publisher = {Elsevier {BV}},
    author = {B. Flemisch and others},
    title = {{DuMux}: {DUNE} for multi-$\lbrace$phase, component, scale, physics, {\ldots}$\rbrace$ flow and transport in porous media},
    journal = {Advances in Water Resources}}
```


# Automated Testing / Test suite

* Latest release (3.8): [![release build badge](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/badges/releases/3.8/pipeline.svg)](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/pipelines?page=1&scope=all&ref=releases/3.8)
* Master branch (development / unstable): [![master build badge](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/badges/master/pipeline.svg)](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/pipelines?page=1&scope=all&ref=master)


DuMu<sup>x</sup> features many tests (some unit tests and test problems) that
are continuously and automatically tested in the GitLab-CI framework (see badges).
The test suite is based on CTest and can also be built and run manually
```
make -j8 build_tests && ctest -j8
```
The tests are labelled (see `CMakeLists.txt` of each individual test for its labels).
You can build and run tests of a specific label (e.g. `2p` for two-phase flow porous medium model tests) like this
```
make -j8 build_2p_tests && ctest -j8 -L ^2p$
```

## Running individual tests

To find out how to build a test inspect the `CMakeLists.txt` file in the respective test folder.

The `dumux_add_test`
command specifies some important parameters: `NAME` sets the name of the test. There is either `SOURCES` or `TARGET`
specified. If `SOURCES` is specified `NAME`, corresponds to the build target, otherwise `TARGET` is the build target.

You can build the test by running `make REPLACE_BY_NAME_OF_BUILD_TARGET` (when using GNU Makefiles (default) this is possible
within the test folder in the build directory, for `ninja` it has to be executed in the top-most build folder level).

Some tests may depend on additional optional dependencies. You can find this by inspecting the argument `CMAKE_GUARD`,
e.g. `HAVE_UMFPACK` means UMFPack is required (via installing Suitesparse), or `( "dune-foamgrid_FOUND" AND "dune-alugrid_FOUND" )`
means that the test requires the additional Dune modules `dune-foamgrid` and `dune-alugrid`. For installing
external dependencies, have a look at the [documentation](https://dumux.org/docs/doxygen/master/external-libraries.html)
and the script [dumux/bin/installexternal.py](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/bin/installexternal.py).

## Test coverage

[![coverage report](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-coverage/badges/master/coverage.svg)](https://pages.iws.uni-stuttgart.de/dumux-repositories/dumux-coverage/)

A weekly coverage report of the test suite is created by gcovr/gcov. The report
currently doesn't include non-instantiated code, so the real coverage is likely lower. However,
only a few lines of code are never instantiated in the comprehensive test suite.


# Contributing

Contributions are highly welcome. Please ask questions over the [mailing list](mailto:dumux@listserv.uni-stuttgart.de).
Please review the [contribution guidelines](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CONTRIBUTING.md)
before opening issues and merge requests.

# Bug/issue reports or vulnerabilities

For bug reports or to report any detected security vulnerabilities contact us
over the mailing list, or file an [issue](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/issues). For bug fixes,
feature implementations open a [merge request](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/merge_requests)
or send us formatted patches.

# Releases and backwards compatibility policy

For a detailed description of the backwards compatibility policy,
please see [contribution guidelines](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CONTRIBUTING.md).

DuMu<sup>x</sup> releases are split into major (e.g. 2.0, 3.0) and minor (e.g. 3.1, 3.2, 3.3) releases.
Major releases are not required to maintain backwards compatibility (see below),
but would provide a detailed guide on how to update dependent modules.
For each minor release, maintaining backwards compatibility is strongly encouraged and recommended.

Despite the goal of maintaining backwards compatibility across minor releases,
for more complicated changes, this is decided upon on a case to case basis, due to limited developer resources.
In the case that implementing full backwards compatibility for an update is not feasible, or would require unreasonable resources,
the degree of backwards compatibility be decided by a vote in one of the monthly core developer meetings.

[0]: https://dumux.org
[1]: https://dune-project.org/
[2]: https://dumux.org/documents/dumux_awrpaper.pdf
[3]: https://dumux.org/docs/doxygen/master/installation.html
[4]: https://www.gnu.org/licenses/gpl-3.0.en.html
[5]: https://www.iws.uni-stuttgart.de/en/lh2/
[6]: https://doi.org/10.1016/j.camwa.2020.02.012
