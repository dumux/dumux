# Part 1: Main program flow

The main program flow is implemented in file `main.cc`, described below. It runs the upscaling
procedure in two stages: `runStaticProblem` drives the quasi-static two-phase invasion-percolation
sweep (`PNMProblemStatic`, see [`problem_static.hh`](../problem_static.hh)) through a full
drainage-imbibition pc cycle, caching the resulting per-throat phase conductances at every step;
`runSinglePhaseUpscaling` then re-solves single-phase Darcy flow once per step and per fluid on
that step's frozen conductance field, converting the resulting mass flux into an apparent
permeability. For the underlying equations and the algorithm in pseudocode, see [Mathematical and
numerical model](../README.md#mathematical-and-numerical-model) in the main description.

The code documentation is structured as follows:

[[_TOC_]]
