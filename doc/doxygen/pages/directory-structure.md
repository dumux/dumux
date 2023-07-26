# Directory Structure

## Top Level

DuMux has the following folder structure, which is similar to other DUNE modules.

* `bin`: executable scripts (Python/Bash) and helper tools, e.g. used for the automatic testing, post-processing, installation
* `cmake`: the configuration options and build system files
* @refdir{doc} "doc": files necessary for the Doxygen documentation and various logos
* @refdir{dumux} "dumux": the main folder, containing the source files. See below for more details.
* [`examples`](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/examples#open_file_folder-example-1-diffusion-equation): well-documented examples of applying DuMux to typical simulation scenarios of different complexity. The example documentation is best viewed with [a browser via GitLab](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/examples#open_file_folder-example-1-diffusion-equation). In the `README.md` file of
each example, the setup is explained and the code is presented and described in detail.
* `test`: unit tests, integration tests and systems tests for each numerical model. The folder structure mimics the structure of the `dumux` folder (see below), the `references` folder contains solutions for the automatic testing. Each system (end user) test describes a full simulation. The setup usually consists of a main file (`main.cc`), the problem definition (`problem.hh`) specifying initial and boundary conditions, and a runtime parameter input file `params.input`. Often spatially-dependent parameters are defined in a separate header (`spatialparams.hh`).

## The Folder dumux

* @refdir{dumux/adaptive} "adaptive": Data transfer between grid views, adaptation indicators.
* @refdir{dumux/assembly} "assembly": Matrix assembler and residual calculation for all discretization methods.
* @refdir{dumux/common} "common": Property system, base classes, boundary conditions, time stepping, splines, dimensionless numbers, type traits, ...
* @refdir{dumux/discretization} "discretization": Infrastructure for discretizations (box, cell-centered, staggered, ...).
* @refdir{dumux/experimental} "experimental": New features, may undergo disruptive changes.
* @refdir{dumux/flux} "flux": Calculation of advective and diffusive fluxes for different discretization schemes.
* @refdir{dumux/freeflow} "freeflow": Single-phase free flow models using Navier-Stokes and eddy-viscosity based Reynolds-averaged Navier-Stokes turbulence models, and shallow water equation model.
* @refdir{dumux/geomechanics} "geomechanics": Elastic and poro-elastic geomechanics models.
* @refdir{dumux/geometry} "geometry": Bounding boxes, intersections, distances, ...
* @refdir{dumux/io} "io": In-/output functionalities such as restart files, gnuplot interface, VTKWriter extensions and grid managers.
* @refdir{dumux/linear} "linear": Linear solver backends.
* @refdir{dumux/material} "material": Constitutive relations and their parameters, definition of components and fluid/solid phases.
* @refdir{dumux/multidomain} "multidomain": Common infrastructure to couple multiple domains of possibly different physics, dimensions and/or locations.
* @refdir{dumux/nonlinear} "nonlinear": Newton's method.
* @refdir{dumux/parallel} "parallel": Helper files for parallel simulations.
* @refdir{dumux/porenetwork} "porenetwork": Models describing a porous medium as a set of pore bodies interconnected by pore throats.
* @refdir{dumux/porousmediumflow} "porousmediumflow": Models for describing flow and mass/momentum/energy transport in a porous medium on the Darcy scale.
* @refdir{dumux/python} "python": Definition of Python bindings for C++ functionalities.
