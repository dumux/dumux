# Benchmarks

The test suite contains a number of verifcation benchmarks
(such as convergence tests against analytical or manufactured solutions).
Moreover, DuMux has been employed in multiple collaborative benchmark projects.
Below, we provide link to descriptions of these benchmarks and instructions
to reproduce the results.

@note This doc page is currently incomplete.
Help exanding the list of documented benchmarks: see [issue tracker](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/issues/1417)
for more information.

## Benchmarks in the test suite

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| ![richards](https://dumux.org/images/mod-richards.png) | @ref benchmark-richards-equation-evaportion "↗️ Richards Evaporation"  | 1D, soil evaporation, @ref RichardsModel, @ref CCTpfaModel | A verfication benchmark with (approximate) analytical solution from Vanderborght et al (2005) @cite Vanderborght2005 |
| ![richards](https://dumux.org/images/mod-richards.png) | @ref benchmark-richards-equation-infiltration "↗️ Richards Infiltration"  | 1D, soil column infiltration front, @ref RichardsModel, @ref CCTpfaModel | A verfication benchmark with (approximate) analytical solution from Vanderborght et al (2005) @cite Vanderborght2005 |
| ![1d rotsym](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/examples/1protationsymmetry/img/result.png) | @ref benchmark-1protsym-example "↗️ Laplace Annulus Rotational Symmetry"  | rotational symmetry, @ref OnePModel | A simple annulus domain with known analytical solution for the Laplace equation with Dirichlet boundary conditions |




## Benchmarks in external modules

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| [![SPE11 Case B CO2](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11/-/raw/main/doc/co2_mass_fraction_liq.gif)](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11) | [↗️ SPE11](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11) | CO<sub>2</sub> Storage, @ref TwoPNCModel, @ref CO2Model  | A DuMux module for running the benchmark scenarios of the [11th Society of Petroleum Engineers' Comparative Solution Project (SPE11) from 2024 (organized by J.M. Nordbotten, M. Fernø, A.R. Kovscek, K.-A. Lie, B. Flemisch)](https://www.spe.org/en/csp/). |
| [![SPE10 Case 2 Fluvial](https://www.spe.org/web/csp/images/img2.gif)](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe10) | [↗️ SPE10](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe10) | Upscaling, @ref TwoPModel | A DuMux module for running the benchmark scenarios of the [10th Society of Petroleum Engineers' Comparative Solution Project (SPE10) from 2001 (organized by M. Christie and M. Blunt)](https://www.spe.org/web/csp/index.html). |
