# Part 1: Implementation of the shallow water flow simulation setup

The shallow water flow setup, including the bottom friction,
is implemented in the files `properties.hh`,
`problem.hh` and `spatialparams.hh`. In the first of these files, a new
type tag is declared for this problem. This allows the specialization
of DuMu<sup>x</sup> `properties` for this type tag, which can be used to customize compile-time
settings for the simulation. Two exemplary `properties`, that are mandatory to be
specialized, are `Problem` and `SpatialParams`. With the first, one sets the
`Problem` class to be used, in which users can define initial and boundary conditions.
Similarly, in the `SpatialParams` class one implements the parameter distributions
(e.g. friction value) that should be used by the model.

The documentation provided in the sequel is structured as follows:

[[_TOC_]]

