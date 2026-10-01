# Benchmarks

The test suite contains a number of verification benchmarks
(such as convergence tests against analytical or manufactured solutions).
Moreover, DuMux has been employed in multiple collaborative benchmark projects.
Below, we provide link to descriptions of these benchmarks and instructions
to reproduce the results.

@note This doc page is currently incomplete.
Help expanding the documentation of implemented benchmarks: see [issue tracker](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/issues/1417)
for more information.

## Porous medium flow benchmarks {#porous-medium-flow-benchmarks}

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| ![richards](https://dumux.org/images/mod-richards.png) | @ref benchmark-richards-benchmarks "↗️ Richards Infiltration (M2.1)" | 1D, soil infiltration, @ref RichardsModel, @ref CCTpfaDiscretization | Travelling-wave benchmark for three soil types (sand, loam, clay); analytical solution from Vanderborght et al. (2005) @cite Vanderborght2005 and Schnepf et al. (2020) @cite Schnepf2020 |
| ![richards](https://dumux.org/images/mod-richards.png) | @ref benchmark-richards-benchmarks "↗️ Richards Evaporation (M2.2)" | 1D, soil evaporation, @ref RichardsModel, @ref CCTpfaDiscretization | Desorptivity-based evaporation benchmark for three soil types; analytical solution from Vanderborght et al. (2005) @cite Vanderborght2005 and Schnepf et al. (2020) @cite Schnepf2020 |
| ![1d rotsym](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/examples/1protationsymmetry/img/result.png) | ↗️ Laplace Annulus Rotational Symmetry (TODO link/readme)  | rotational symmetry, @ref OnePModel | A simple annulus domain with known analytical solution for the Laplace equation with Dirichlet boundary conditions |
| ![buckley-leverett](buckleyleverett_lineplot_comparison.png) | @ref benchmark-buckley-leverett "↗️ Buckley-Leverett" | @ref TwoPModel, @ref CCTpfaDiscretization | Two-phase fluid displacement |
| ![mcwhorter-sunada](mcwhortersunada_lineplot_comparison.png) | @ref benchmark-mcwhorter-sunada "↗️ McWhorter-Sunada" | @ref TwoPModel, @ref CCTpfaDiscretization | Counter-current imbibition compared with Fučík's semi-analytical reference |
| ![horton-rogers-lapwood](hortonrogerslapwood_cells.png) | @ref benchmark-horton-rogers-lapwood "↗️ Onset of Convection (Horton-Rogers-Lapwood)" | density-driven flow, Boussinesq approximation, @ref OnePNCModel, @ref BoxDiscretization | Growth rate of a convection cell in a porous layer with dense fluid on top, compared with linear stability theory @cite HortonRogers1945 @cite Lapwood1948 |
| ![heat-pipe](heatpipe2p1c_profile.png) | @ref benchmark-heatpipe-udell "↗️ Udell's Heat Pipe" | 1D, non-isothermal, phase change, @ref TwoPOneCModel, @ref BoxDiscretization, @ref CCTpfaDiscretization | Steady counter-current liquid-vapor flow driven by a heat flux; first-order convergence with full upwinding and second-order convergence with central mobility weighting towards semi-analytical solution of Udell (1985) @cite Udell1985HeatPipe |
| ![1pni convection](1pni_1d_convection_benchmark_lineplot.png) | @ref benchmark-1pni-convection "↗️ 1D Heat Convection with a Retarded Thermal Front" | 1D, advective heat transport, @ref OnePModel, @ref NIModel, @ref CCTpfaDiscretization, @ref CCMpfaDiscretization, @ref BoxDiscretization | Retarded thermal front in a fully saturated porous medium; front position from the barycenter against the retarded step front, front shape by grid convergence against the advection-diffusion solution of Ogata and Banks (1961) @cite OgataBanks1961 |
| [![SPE11 Case B CO2](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11/-/raw/main/doc/co2_mass_fraction_liq.gif)](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11) | [↗️ SPE11](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe11) | CO<sub>2</sub> Storage, @ref TwoPNCModel, @ref CO2Model  | A DuMux module for running the benchmark scenarios of the [11th Society of Petroleum Engineers' Comparative Solution Project (SPE11) from 2024 (organized by J.M. Nordbotten, M. Fernø, A.R. Kovscek, K.-A. Lie, B. Flemisch)](https://www.spe.org/en/csp/). |
| [![SPE10 Case 2 Fluvial](https://www.spe.org/web/csp/images/img2.gif)](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe10) | [↗️ SPE10](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-spe10) | Upscaling, @ref TwoPModel | A DuMux module for running the benchmark scenarios of the [10th Society of Petroleum Engineers' Comparative Solution Project (SPE10) from 2001 (organized by M. Christie and M. Blunt)](https://www.spe.org/web/csp/index.html). |
| ![radial injection](2pve_radialinjection_profile.png) | @ref benchmark-2pve-radial-injection "↗️ Radial Injection into a Confined Aquifer" | CO<sub>2</sub> storage, vertical equilibrium, @ref TwoPVEModel, @ref CCTpfaDiscretization | Interface of a CO<sub>2</sub> plume injected into a confined aquifer; comparison with the similarity solution of Nordbotten and Celia (2006) @cite NordbottenCelia2006 |

## Free flow benchmarks {#free-flow-benchmarks}

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| [![lid-driven cavity](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/examples/liddrivencavity/img/result.svg)](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/liddrivencavity/README.md) | [↗️ Shear-Driven Cavity Flow](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/liddrivencavity/README.md) | 2D, free flow, @ref NavierStokesModel, @ref FaceCenteredStaggeredDiscretization | Velocity profiles at Re = 1 and Re = 1000 compared against reference data from Ghia et al. (1982) and Jurjević (1999) |
| ![taylor-couette](taylorcouette_comparison.png) | @ref benchmark-taylor-couette "↗️ Taylor-Couette Flow" | 2D, free flow, @ref NavierStokesModel, @ref DiamondDiscretization, @ref CCTpfaDiscretization | Circular Couette flow between two rotating cylinders at Re = 100; second-order convergence of velocity and pressure towards the analytical solution @cite taylor1923 |

## Shallow water benchmarks {#shallow-water-benchmarks}

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| ![dam break](shallowwater_dambreak_setup.svg) | @ref benchmark-dam-break "↗️ Wet Dam Break" | Riemann problem, @ref ShallowWaterModels, @ref CCTpfaDiscretization | Bore and rarefaction wave after the removal of a gate; the relative $L^1$ errors converge with approximately first order towards the exact solution @cite Stoker1957, with the same rates as a textbook first-order Godunov scheme |
| ![paraboloid bowl](shallowwater_bowl_setup.svg) | @ref benchmark-thacker-bowl "↗️ Oscillation in a Paraboloid Bowl (Thacker)" | wetting and drying, @ref ShallowWaterModels, @ref CCTpfaDiscretization | Oscillation of the free surface in a bowl without friction with a moving shoreline; the first-order scheme damps the oscillation, and the relative $L^1$ error of the water depth converges with approximately first order on fine grids towards the exact solution @cite Thacker1981 |
| ![lake at rest](shallowwater_lakeatrest_setup.svg) | @ref benchmark-lake-at-rest "↗️ Lake at Rest" | well-balancing, wetting and drying, @ref ShallowWaterModels, @ref CCTpfaDiscretization | Water at rest over a hump @cite LeVeque1998 and a sill with discontinuous flanks; preserved to round-off on a submerged bed |
| ![Poiseuille flow](shallowwater_poiseuille_setup.svg) | @ref benchmark-shallow-water-poiseuille "↗️ Shallow Water Poiseuille Flow" | viscous flux, bed friction, @ref ShallowWaterModels, @ref CCTpfaDiscretization | Steady channel flow with an eddy viscosity and no-slip side walls, and laminar film flow with viscous bed friction; approximately first-order convergence along the channel, caused by the bed being constant in each cell, and approximately second-order convergence across the channel towards the exact velocity |
| ![rough channel](shallowwater_roughchannel_setup.svg) | @ref benchmark-rough-channel "↗️ Flow in a Rough Channel" | bed friction, @ref ShallowWaterModels, @ref CCTpfaDiscretization | Normal flow with the friction law of Gauckler, Manning and Strickler; approximately first-order convergence towards the normal depth, caused by the bed being constant in each cell |

## Solid mechanics benchmarks {#solid-mechanics-benchmarks}

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| ![membrane plate](plate_membrane_deformation_3d.png) | @ref benchmark-membrane-plate "↗️ Clamped Circular Membrane" | plate mechanics, @ref MembranePlate, @ref BoxDiscretization | Clamped circular membrane under uniform load; convergence against analytical solution |
| ![kirchhoff-love plate](plate_kirchhoff_love_deformation_3d.png) | @ref benchmark-kirchhoff-love-plate "↗️ Clamped Circular Plate (Kirchhoff-Love)" | plate mechanics, @ref KirchhoffLovePlate, @ref BoxDiscretization, @ref PQ1BubbleDiscretization, @ref MultiDomain | Clamped circular plate under uniform load; convergence against analytical solution |
| ![mindlin-reissner plate](plate_mindlin_reissner_deformation_3d.png) | @ref benchmark-mindlin-reissner-plate "↗️ Clamped Circular Plate (Mindlin-Reissner)" | plate mechanics, @ref MindlinReissnerPlate, @ref BoxDiscretization, @ref PQ1BubbleDiscretization, @ref MultiDomain | Clamped circular plate under uniform load with shear correction; convergence against analytical solution |
|  | @ref benchmark-cooks-membrane-elastic "↗️ Cook's Membrane, Linear Elasticity" | solid mechanics, linear elasticity, locking, @ref Elastic, @ref BoxDiscretization, PQ1-bubble and PQ2 | Tapered panel under shear load with a nearly incompressible material; tip displacement under mesh refinement compared with @cite Chen2024 |
| ![hyperelastic](solid_turekhron_domain.svg) | @ref benchmark-turekhron-csm "↗️ Hyperelastic Dynamic (CSM1-3)" | solid mechanics, hyperelasticity, @ref BoxDiscretization | oscillating/flapping thin solid rectangle ; comparison against reference data from @cite TurekHron2006 |
| ![cook hyperelastic](cooks_membrane_hyperelastic_deformation.png) | @ref benchmark-cooks-membrane-hyperelastic "↗️ Cook's Membrane, Hyperelastic (SPP 1748)" | solid mechanics, hyperelasticity, locking, mixed formulation, @ref HyperelasticVolIso, MINI and Taylor–Hood elements | Tapered cantilever under shear load, 2D plane strain and 3D, Poisson's ratio from 0.35 to 0.4999; tip displacement compared with @cite Schroeder2021 |
| ![incompressible block](incompressible_block_deformation.png) | @ref benchmark-incompressible-block "↗️ Incompressible Block (SPP 1748)" | solid mechanics, hyperelasticity, locking, mixed formulation, @ref HyperelasticVolIso, MINI, Taylor–Hood and Q2–Q1 elements | Nearly incompressible block under a partial surface load; displacement of the load centre compared with @cite Schroeder2021 |
| [![cardiac elastodynamics](https://bwgit.de/st/iws/dumux/appl/dumux-cardiac-benchmark/-/raw/main/doc/media/biventricle_heartbeat_thumbnail.gif)](https://st.bwpages.de/iws/dumux/appl/dumux-cardiac-benchmark/) | [↗️ Cardiac Elastodynamics](https://st.bwpages.de/iws/dumux/appl/dumux-cardiac-benchmark/) | solid mechanics, hyperelasticity, elastodynamics, @ref Hyperelastic, PQ2 finite elements | A DuMux module for the heartbeat of an idealized left ventricle and biventricle (anisotropic, nearly incompressible, viscous, actively contracting material); probe displacements compared with the participants of @cite Arostica2025 ([code](https://bwgit.de/st/iws/dumux/appl/dumux-cardiac-benchmark)) |

## Poromechanics benchmarks {#poromechanics-benchmarks}

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| ![mandel](mandel_pressure.png) | @ref benchmark-mandel "↗️ Mandel's Problem" | poromechanics, consolidation, @ref OnePModel, @ref PoroElastic, @ref CCTpfaDiscretization, @ref BoxDiscretization, @ref MultiDomain | Consolidation of a poroelastic slab between rigid plates (Mandel-Cryer effect), also with compressible solid grains; comparison with the analytical solution @cite Abousleiman1996 |
| ![cryer sphere](cryer_overview.png) | @ref benchmark-cryer "↗️ Cryer Sphere Consolidation" | poromechanics, consolidation, large deformation, @ref OnePModel, @ref PoroElastic, @ref PoroElasticLargeDef, @ref CCMpfaDiscretization, @ref PQ2Discretization, @ref BoxDiscretization, @ref MultiDomain | Mandel–Cryer effect: the centre pore pressure of a drained poroelastic sphere overshoots the surface load; linear and large-deformation poroelasticity, comparison with Cryer's solution @cite Cryer1963 |

## Time integration benchmarks {#time-integration-benchmarks}

| Image | Benchmark description | Topics | Comment |
|:----------:|:----------:|:---------|:---------|
| ![time-stepping stability](timestepping_stability_comparison.png) | @ref benchmark-timestepping-methods "↗️ Time-Stepping Methods" | time integration, `MultiStageMethod`, `MultiStageTimeStepper` | Observed temporal convergence order and linear stability regions of 8 explicit, implicit and DIRK schemes, verified against analytical references |
