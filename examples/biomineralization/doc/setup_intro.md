# Simulation Setup

The setup for this example is based on the initial 100000s of the _column experiment D1_ described in [@Hommel2015] and can also be found in Section 4.2 of [@Hommel2016] under the name _column experiment D1_.
We here only consider the "simplified chemistry case" in which it is assumed that the precipitation rate is equal to the ureolysis rate to simplify the chemical processes considered ($`r_\mathrm{prec} = r_\mathrm{urea}`$) (see Chapter 6 [@Hommel2016]).

The column is considered as a 1D domain with injections of various aqueous solutions at the bottom (Neumann boundary condition)
and a fixed pressure at the top (Dirichlet boundary condition).
The various types of injected fluid compositions are read from an additional input file `injections_type.dat`.
Depending on the injection type, the `problem.hh` adapts the composition of the injected aqueous solution for the bottom Neumann boundary condition.
Further, `spatialparams.hh` is an example for spatial parameters dealing with changing porosity and permeability.

The subsequent file documentation is structured as follows:

[[_TOC_]]

[@Hommel2015]: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014WR016503 "A revised model for microbially induced calcite precipitation: Improvements	and new insights based on recent experiments"
[@Hommel2016]: https://elib.uni-stuttgart.de/handle/11682/8787 "Modelling biogeochemical and mass transport processes in the subsurface: investigation of microbially induced calcite precipitation"
