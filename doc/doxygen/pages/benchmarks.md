# Benchmarks

The test suite contains a number of verification benchmarks
(such as convergence tests against analytical or manufactured solutions).
Moreover, DuMux has been employed in multiple collaborative benchmark projects.
Below, we provide link to descriptions of these benchmarks and instructions
to reproduce the results.

@note This doc page is currently incomplete.
Help expanding the documentation of implemented benchmarks: see [issue tracker](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/issues/1417)
for more information.

## Benchmarks in the test suite

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| ![richards](https://dumux.org/images/mod-richards.png) | @ref benchmark-richards-benchmarks "↗️ Richards Infiltration (M2.1)" | 1D, soil infiltration, @ref RichardsModel, @ref CCTpfaDiscretization | Travelling-wave benchmark for three soil types (sand, loam, clay); analytical solution from Vanderborght et al. (2005) @cite Vanderborght2005 and Schnepf et al. (2020) @cite Schnepf2020 |
| ![richards](https://dumux.org/images/mod-richards.png) | @ref benchmark-richards-benchmarks "↗️ Richards Evaporation (M2.2)" | 1D, soil evaporation, @ref RichardsModel, @ref CCTpfaDiscretization | Desorptivity-based evaporation benchmark for three soil types; analytical solution from Vanderborght et al. (2005) @cite Vanderborght2005 and Schnepf et al. (2020) @cite Schnepf2020 |
| ![1d rotsym](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/examples/1protationsymmetry/img/result.png) | ↗️ Laplace Annulus Rotational Symmetry (TODO link/readme)  | rotational symmetry, @ref OnePModel | A simple annulus domain with known analytical solution for the Laplace equation with Dirichlet boundary conditions |
| ![membrane plate](plate_membrane_deformation_3d.png) | @ref benchmark-membrane-plate "↗️ Clamped Circular Membrane" | plate mechanics, @ref MembranePlate, @ref BoxDiscretization | Clamped circular membrane under uniform load; convergence against analytical solution |
| ![kirchhoff-love plate](plate_kirchhoff_love_deformation_3d.png) | @ref benchmark-kirchhoff-love-plate "↗️ Clamped Circular Plate (Kirchhoff-Love)" | plate mechanics, @ref KirchhoffLovePlate, @ref BoxDiscretization, @ref PQ1BubbleDiscretization, @ref MultiDomain | Clamped circular plate under uniform load; convergence against analytical solution |
| ![mindlin-reissner plate](plate_mindlin_reissner_deformation_3d.png) | @ref benchmark-mindlin-reissner-plate "↗️ Clamped Circular Plate (Mindlin-Reissner)" | plate mechanics, @ref MindlinReissnerPlate, @ref BoxDiscretization, @ref PQ1BubbleDiscretization, @ref MultiDomain | Clamped circular plate under uniform load with shear correction; convergence against analytical solution |
| ![buckley-leverett](buckleyleverett_lineplot.png) | @ref benchmark-buckley-leverett "↗️ Buckley-Leverett" | @ref TwoPModel, @ref CCTpfaDiscretization | Two-phase fluid displacement |

## Benchmarks in documented examples

| Image |                                Benchmark description                                 | Topics | Comment |
|:------------------------------------------------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------:|:---------|:---------|
| [![lid-driven cavity](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/examples/liddrivencavity/img/result.svg)](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/liddrivencavity/README.md) | [↗️ Shear-Driven Cavity Flow](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/liddrivencavity/README.md) | 2D, free flow, @ref NavierStokesModel, @ref FaceCenteredStaggeredDiscretization | Velocity profiles at Re = 1 and Re = 1000 compared against reference data from Ghia et al. (1982) and Jurjević (1999) |

## Benchmarks in external modules

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| [![SPE11 Case B CO2](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11/-/raw/main/doc/co2_mass_fraction_liq.gif)](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11) | [↗️ SPE11](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11) | CO<sub>2</sub> Storage, @ref TwoPNCModel, @ref CO2Model  | A DuMux module for running the benchmark scenarios of the [11th Society of Petroleum Engineers' Comparative Solution Project (SPE11) from 2024 (organized by J.M. Nordbotten, M. Fernø, A.R. Kovscek, K.-A. Lie, B. Flemisch)](https://www.spe.org/en/csp/). |
| [![SPE10 Case 2 Fluvial](https://www.spe.org/web/csp/images/img2.gif)](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe10) | [↗️ SPE10](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe10) | Upscaling, @ref TwoPModel | A DuMux module for running the benchmark scenarios of the [10th Society of Petroleum Engineers' Comparative Solution Project (SPE10) from 2001 (organized by M. Christie and M. Blunt)](https://www.spe.org/web/csp/index.html). |
