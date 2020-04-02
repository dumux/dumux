# Part 1: Implementation of the single-phase flow simulation setup

The single-phase flow setup is implemented in the files `properties_1p.hh`,
`problem_1p.hh` and `spatialparams_1p.hh`. In the first of these files, a new
type tag is declared for this problem. This then allows the specialization
of DuMuX `properties` for this type tag, which can be used to customize compile-time
settings for the simulation. Two exemplary `properties`, that are mandatory to be
specialized, are `Problem` and `SpatialParams`. With the first, one sets the
`Problem` class to be used, in which users can define initial and boundary conditions.
Similarly, in the `SpatialParams` class one implements the parameter distributions
(e.g. porosity and permeability) that should be used by the model.

The documentation provided in the sequel is structured as follows:

[[_TOC_]]
