# Main file

The main file features a time loo with check points to assure the correct timing of the injections matching an experimental setup (see [@Hommel2016], Section 4.2 under the name _column experiment D1_)
The timing of the injections is stored in an additional input file `injections_checkpoints.dat`,
which is read by `main.cc` to set the check points in the time loop.
The number of check points reached is counted and passed on to `problem.hh`, so that `problem.hh` can determine the appropriate injection type.

Additionally, there is an output of the porosity-dependent permeability.

The subsequent file documentation is structured as follows:

[[_TOC_]]

[@Hommel2016]: https://elib.uni-stuttgart.de/handle/11682/8787 "Modelling biogeochemical and mass transport processes in the subsurface: investigation of microbially induced calcite precipitation"
