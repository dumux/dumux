# Examples and tutorials

You have downloaded DuMux and its dependencies.
You have run `dunecontrol` and your first example compiled and showed a nice simulation in ParaView.
What now? *How on earth is this going to help me solve my multi-(phase, component, scale, physics) flow and transport problems in porous media systems?*
A great collection of additional resources can be found below.

| | | |
|:----------------:|:---------:|:--------------|
| <div class="dumux-circled-text">📝</div>   | [Documented Examples](@ref examples)  | DuMux comes with a set of [➡️ Documented Examples](@ref examples). These are currently best viewed in GitLab via the browser. The corresponding files are located in the same folder as the `README.md` which contains the documentation and is rendered as a documentation page in GitLab. |
| <div class="dumux-circled-text">🧪</div>   | [↗️ DuMux test applications](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/test)                  |  Every model in DuMux is accompanied by at least one test application. Applications can be used as a starting point for developing your own applications. All tests can be found in the DuMux source code folder [↗️ dumux/test](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/test) and instruction how to build and run tests is described [here](#running-individual-tests). |
| <div class="dumux-circled-text">🧑‍🏫</div>   | [↗️ DuMux course](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course)                  | This dumux course is a 3 day workshop offered in person from time to time. But all the course materials is also [↗️ available online for self study](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course). A series of beginner-level exercises are explained such that you can see how a model is developed in DuMux. There is exercises and solution on models exploring "Coupling free flow and porous-media flow", "Flow in fractured porous media" and "Fluid-solid phase change". Installation instructions and further instructions how to use the course material can be found [↗️ here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course). |
| <div class="dumux-circled-text">🎓</div>   | [↗️ DuMux lecture](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-lecture)                  | This repository contains different application examples that are used in the university-level teaching at the Department of Hydromechanics and Modelling of Hydrosystems at the University of Stuttgart. See more information below. |


__Details on the dumux lecture repository__

The `dumux-lecture` module can be
installed via the [`bin/installexternal.py` script](#install-external-dependencies-via-script).

The module is structured based on the different lectures:

* `mm`: Multiphase Modelling
* `efm`: Environmental Fluid Mechanics
* `mhs`: Modelling of Hydrosystems

The majority of applications are covered in the course "Multiphase Modelling" (`mm` folder),
while there are also some basic examples in the courses "Environmental Fluid Mechanics" (`efm` folder)
and "Modelling of Hydrosystems" (`mhs` folder).

These applications are primarily designed to enhance the understanding of conceptualizing the governing physical processes and their implementation in a numerical simulator. Different aspects of modeling multi-phase multi-component flow and transport processes are shown. The lectures focus on questions such as the assignment of boundary conditions, the choice of the appropriate physics for a given problem (which phases, which components), discretization issues, time stepping. You can find, for example, a comparison of different two-phase flow problems: The simpler approach considers two immiscible fluids while components in both phases with inter-phase mass transfer are considered in the more complex approach. All scenarios and their physical background are explained in additional LaTeX (`*.tex`) files, which are provided in sub-directories named description. The following test cases are contained in the `dumux-lecture` module:

* `buckleyleverett`: The Buckley-Leverett Problem is a classical porous media flow show case
* `co2plume`: Analysis of the influence of the gravitational number on a CO2 plume
* `columnxylene`: An experiment of the Research Facility for Subsurface Remediation, University of Stuttgart
* `convectivemixing`: A test case related to CO2 storage
* `fractures`: Two-phase flow in fractured porous media
* `fuelcell`: Water management in PEM fuel cells
* `heatpipe`: A show case for two-phase two-component flow with heat fluxes
* `heavyoil`: Steam assisted gravity drainage (SAGD)
* `henryproblem`: A show case related to salt water intrusion
* `mcwhorter`: The McWhorter Problem is a classical porous media flow show case
* `naplinfiltration`: Infiltration of non-aqueous phase liquid (NAPL) into soil
* `remediationscenarios`: Test case for NAPL contaminated unsaturated soils
* `groundwater`: Simple groundwater flow case for the course Modelling of Hydrosystems (mhs)
