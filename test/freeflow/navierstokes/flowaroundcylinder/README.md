# DFG Flow Around a Cylinder (Schäfer & Turek) {#benchmark-dfg-flow-around-cylinder}

This test implements the laminar flow around a cylinder benchmark of the DFG priority research
program "Flow Simulation with High-Performance Computers" @cite SchaeferTurek1996.

The benchmark defines a two-dimensional channel flow around a circular obstacle. Highly accurate
reference values are available for three quantities: the drag coefficient, the lift coefficient
and a pressure difference.

All three 2D cases are implemented:

* the stationary case 2D-1 ($\text{Re} = 20$) for all discretization variants;
* the time-dependent cases 2D-2 and 2D-3 ($\text{Re} = 100$) for the three PQ1Bubble/Box variants.

## Problem description

The incompressible Navier-Stokes equations

$$
\begin{aligned}
\rho\left(\frac{\partial \mathbf{v}}{\partial t} + \nabla\cdot(\mathbf{v}\otimes\mathbf{v})\right)
  - \mu\Delta\mathbf{v} + \nabla p &= \mathbf{0}, \\
\nabla\cdot\mathbf{v} &= 0
\end{aligned}
$$

are solved in the channel $\Omega = [0, 2.2] \times [0, 0.41] \setminus B_r(0.2, 0.2)$.

The cylinder has the radius $r = 0.05$ (diameter $D = 0.1$). It is placed slightly below the
channel centre line, so the setup is not symmetric and the lift coefficient is non-zero.

![Domain and boundary conditions of the DFG benchmark (to scale).](dfg_cylinder_domain.svg)

| Parameter | Symbol | Value | Unit |
| :--- | :---: | :---: | :---: |
| Density | $\rho$ | $1$ | $\text{kg/m}^3$ |
| Kinematic viscosity | $\nu = \mu/\rho$ | $10^{-3}$ | $\text{m}^2\text{/s}$ |
| Cylinder diameter | $D$ | $0.1$ | $\text{m}$ |
| Channel height | $H$ | $0.41$ | $\text{m}$ |
| Pressure evaluation point (cylinder front) | $A$ | $(0.15, 0.2)$ | $\text{m}$ |
| Pressure evaluation point (cylinder back) | $B$ | $(0.25, 0.2)$ | $\text{m}$ |

### The three 2D cases

The mean inflow velocity $\bar{U} = \tfrac{2}{3}U$ and the cylinder diameter $D$ define the Reynolds number
$\text{Re} = \bar{U}D/\nu$.

| Case | Inflow peak velocity $U(t)$ [m/s] | $\bar{U}$ [m/s] | Re | Character | Parameter file |
| :--- | :---: | :---: | :---: | :--- | :--- |
| 2D-1 | $0.3$ | $0.2$ | $20$ | stationary | `params_1.input` |
| 2D-2 | $1.5$ | $1.0$ | $100$ | periodic (Kármán vortex street) | `params_2.input` |
| 2D-3 | $1.5\sin(\pi t/8)$, $t \in [0, 8]\,\text{s}$ | $1.0$ | $0 \to 100 \to 0$ | time-dependent inflow | `params_3.input` |

The parameters `Problem.InflowMaxVelocity` ($U$), `Problem.InflowTimeDependent` (multiplies $U$ by
$\sin(\pi t/8)$) and `Problem.Instationary` (time loop) select the case.

### Boundary and initial conditions

* **Inflow** ($x = 0$): parabolic velocity profile $\mathbf{v}(0, y, t) = \left(4U(t)\, y (H - y)/H^2,\ 0\right)^T$.
* **Walls** ($y = 0$, $y = H$) and **cylinder**: no-slip, $\mathbf{v} = \mathbf{0}$.
* **Outflow** ($x = 2.2$): do-nothing condition.
  * The advective momentum flux $\rho(\mathbf{v}\cdot\mathbf{n})\mathbf{v}$ leaves the domain.
  * The viscous and pressure contributions of the normal stress are zero.
* **Initial condition** (cases 2D-2 and 2D-3): $\mathbf{v} = \mathbf{0}$, $p = 0$.

## Quantities of interest (benchmark indicators)

The drag and lift forces act on the cylinder surface $S$. Let $\mathbf{n}$ be the unit normal on $S$
pointing into the cylinder and $v_\tau$ the tangential velocity. Then

$$
F_D = \int_S \left(\rho\nu\frac{\partial v_\tau}{\partial \mathbf{n}} n_y - p\, n_x\right) \text{d}S, \qquad
F_L = -\int_S \left(\rho\nu\frac{\partial v_\tau}{\partial \mathbf{n}} n_x + p\, n_y\right) \text{d}S .
$$

The benchmark indicators are the drag coefficient $C_D = 2 F_D/(\rho\bar{U}^2 D)$, the lift coefficient
$C_L = 2 F_L/(\rho\bar{U}^2 D)$ and the pressure difference $\Delta p = p(A) - p(B)$.
For case 2D-3, $\bar{U} = 1$ is the mean velocity of the peak inflow profile (not of $U(t)$).

For case 2D-2, $f$ denotes the frequency of the vortex shedding.

The FeatFlow reference values are listed below (see the references at the end). For 2D-2 and 2D-3, they are
computed with $Q_2/P_1^\text{disc}$ elements on level 6 and the Crank-Nicolson method with $\Delta t = 1/1600\,\text{s}$.

| Case | Quantity | Reference value |
| :--- | :--- | :--- |
| 2D-1 | $C_D$ | $5.57953523384$ |
| 2D-1 | $C_L$ | $0.010618948146$ |
| 2D-1 | $\Delta p$ | $0.11752016697$ |
| 2D-2 | $C_{D,\max}$ | $3.22739$ |
| 2D-2 | $C_{L,\max}$ | $0.98657$ |
| 2D-3 | $C_{D,\max}$ at $t$ | $2.9437637214$ at $t = 3.9365625\,\text{s}$ |
| 2D-3 | $C_{L,\max}$ at $t$ | $0.47748781595$ at $t = 5.6928125\,\text{s}$ |
| 2D-3 | $\Delta p(t = 8\,\text{s})$ | $-0.11154138872$ |

### Evaluation in the code

* **Drag and lift**: the traction is integrated over the boundary faces on the cylinder.
  On the no-slip cylinder wall, this is equivalent to the formula above.
  The viscous stress is the one of the discretization: $\mu \nabla \mathbf{v}$ with
  `FreeFlow.EnableUnsymmetrizedVelocityGradient = true` (set in all input files), otherwise
  $\mu (\nabla \mathbf{v} + \nabla \mathbf{v}^T)$. The input files use the unsymmetrized gradient, since the
  natural (do-nothing) outflow condition of the symmetrized form, $\boldsymbol{\sigma}\mathbf{n} = 0$, changes the
  flow near the outlet (Diamond/CCTpfa: $C_D = 4.72$ instead of $5.57$ on the default grid).
* **Pressure difference**: the pressure is interpolated at $A$ and $B$.
  * The cylinder surface is approximated by a polygon, so $A$ and $B$ may lie outside the grid.
  * In that case, they are moved upstream ($A$) or downstream ($B$) in steps of $1\,\text{mm}$ until
    they hit an element.
* The evaluation is not parallelized. A parallel run with inertia terms stops with an exception.

## Discretization

The momentum balance and the mass balance are discretized separately and coupled by the
`FreeFlowCouplingManager`. Five combinations of discretization schemes and grids are tested:

| Momentum | Mass | Grid |
| :--- | :--- | :--- |
| @ref DiamondDiscretization | @ref CCTpfaDiscretization | triangles |
| @ref PQ1BubbleDiscretization | @ref BoxDiscretization | triangles |
| @ref PQ1BubbleDiscretization (hybrid) | @ref BoxDiscretization | triangles |
| @ref PQ1BubbleDiscretization (hybrid) | @ref BoxDiscretization | quadrilaterals |
| @ref PQ1BubbleDiscretization | @ref DiamondDiscretization | quadrilaterals |

The hybrid PQ1Bubble scheme (`PQ1BubbleHybridModel`) uses the same basis functions. Only the vertex dofs have
control volumes, the bubble dofs are treated with the finite element method. On quadrilaterals, it uses two bubble
functions per element. The quadrilateral variant runs the same executable on the quadrilateral mesh, so it allows
comparing the grids with the same discretization and the same problem interface.

The three PQ1Bubble/Box variants use the new CVFE problem and variables interfaces (`problem_newinterface.hh`,
compile definition `NEW_PROBLEM_INTERFACE=1`). The other two variants use the established interface (`problem.hh`).
The two problem files differ as follows:

| | `problem.hh` | `problem_newinterface.hh` |
| :--- | :--- | :--- |
| Dirichlet conditions | boundary types per control volume | constraints on the boundary dofs |
| Drag and lift | traction per boundary face, gradient integrated with a quadrature rule of order 5 and divided by the face area, pressure at the face integration point | traction evaluated and integrated with a quadrature rule of order 5 |

The advective term is discretized with central differences (`Flux.UpwindWeight = 0.5`).

The unstructured grids are generated with Gmsh. Both use a mesh size of $h = 0.01$ in the channel and
$h = 0.002$ on the cylinder surface.

| Mesh | Elements |
| :--- | :--- |
| `meshes/cylinder_channel.geo` | 28948 triangles |
| `meshes/cylinder_channel_quad.geo` | 13826 quadrilaterals |

## Time discretization (cases 2D-2 and 2D-3)

The time-dependent cases use the second-order DIRK scheme of Alexander @cite Alexander1977 (`TimeLoop.Scheme = DIRK2`, default).
`DIRK3`, `ImplicitEuler` and `CrankNicolson` can be chosen as well.

* DIRK2 is second-order accurate, L-stable and stiffly accurate. First-order schemes such as the implicit Euler
  method damp the oscillations @cite SchaeferTurek1996.
* It is run with the multi-stage assembler and time stepper (`Experimental::MultiStageMultiDomainAssembler`,
  `Experimental::MultiStageTimeStepper`). These are only available for the CVFE schemes with the new problem
  interface, i.e. the three PQ1Bubble/Box variants.
* The time-dependent inflow (Dirichlet constraints) is updated at each stage time.
* Newton's method uses the shift criterion (`Newton.MaxRelativeShift`).
* The Jacobian matrix is reused (`Newton.EnableJacobianReuse = true`) while the Newton updates contract. With DIRK2 and a
  fixed time step size, the stage weights do not change, so a Jacobian matrix can be reused over stages and time steps.
* The linear systems are solved with MUMPS (`DirectSolverMumps`) if available, otherwise with UMFPack.

Fixed time step sizes are used, following the benchmark description:

| Case | Time steps | Parameters |
| :--- | :--- | :--- |
| 2D-2 | coarse step $0.1\,\text{s}$ until $3.5\,\text{s}$ to develop the periodic flow, then $0.001\,\text{s}$ until $10\,\text{s}$ | `TimeLoop.Dt1`, `T1`, `Dt2`, `T2` |
| 2D-3 | $1/1600\,\text{s}$ on $[0, 8]\,\text{s}$ | `TimeLoop.Dt`, `TEnd` |

If Newton's method fails, the step is repeated with a reduced step size. Afterwards, the fixed step size is used again.

### Evaluation at the time levels

The benchmark indicators are evaluated from the solution at the new time level $t^{n+1}$.

* The DIRK schemes do not use the spatial terms of the old time level. The pressure, the Lagrange multiplier of the
  incompressibility constraint, is therefore determined at each time level.
* The Crank-Nicolson method ($\theta = 0.5$) weights all spatial terms with $\tfrac{1}{2}$, including the
  pressure gradient and the continuity equation. Only $\tfrac{1}{2}(p^n + p^{n+1})$ is determined, and the
  pressure at the time levels shows an undamped odd-even oscillation (e.g. $\Delta p$ alternating between
  $-1.9$ and $6.9$ in case 2D-2). With `CrankNicolson`, the pressure-dependent indicators are not meaningful.

The FeatFlow time series are given at the half steps of their Crank-Nicolson method, so the times of their extrema
are shifted by $\Delta t/2$ with respect to the DuMux values.

Output of the time-dependent cases (file names from `Problem.IndicatorFile`, default `Problem.Name`):

| File | Content | Written |
| :--- | :--- | :--- |
| `<name>.csv` | `time`, `cDrag`, `cLift`, `pDiff` of every time step | appended and flushed at every step, so it survives an aborted run (e.g. the wall-time limit of a cluster) |
| `<name>.json` | the same time series plus grid size, time stepping scheme, `finished` and `error` | at the end of the simulation and if it fails with an exception |
| `<Problem.Name>-*.vtu`, `.pvd` | velocity and pressure fields (binary) | every `Vtk.OutputTimeInterval` = 0.01 s, i.e. about 33 frames per shedding period ($1/f \approx 1/3\,\text{s}$), independent of the time step size |

## Running the tests

```bash
cd build-cmake/test/freeflow/navierstokes/flowaroundcylinder
make test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box
./test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params_1.input  # case 2D-1
./test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params_2.input  # case 2D-2 (days)
./test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params_3.input  # case 2D-3 (days)
```

The executable is named after its regression test (case 2D-1), but runs all three cases. The hybrid variant
(`test_ff_navierstokes_dfg_benchmark_stationary_hybrid_pq1bubble_box`) is run in the same way, on quadrilaterals
with `-Grid.File meshes/cylinder_channel_quad.msh` (tests `*_hybrid_pq1bubble_box_quad`).

The ctest regression tests only cover the stationary case 2D-1. The time-dependent cases are run
with `benchmark.py` (see below).

The tests come in two variants:

* `test_ff_navierstokes_dfg_benchmark_stationary_*` solve case 2D-1 and print $C_D$, $C_L$ and
  $\Delta p$ together with the reference values.
* `test_ff_stokes_dfg_benchmark_stationary_*` solve the same setup without the inertia terms
  (`Problem.EnableInertiaTerms = false`). They only serve as regression tests.

`Problem.CheckIndicators = true` makes the run fail if the indicators deviate too much from the reference
values. This only passes on refined grids.

## Results on the default grids

| Momentum / mass | $C_D$ | $C_L$ | $\Delta p$ |
| :--- | :---: | :---: | :---: |
| Diamond / CCTpfa | $5.5664$ | $0.00949407$ | $0.115135$ |
| PQ1Bubble / Box | $5.51483$ | $0.00896245$ | $0.118022$ |
| PQ1Bubble (hybrid) / Box | $5.5236$ | $0.00892858$ | $0.118087$ |
| PQ1Bubble (hybrid) / Box (quadrilaterals) | $5.53695$ | $0.00669747$ | $0.117841$ |
| PQ1Bubble / Diamond (quadrilaterals) | $5.61234$ | $0.00387379$ | $0.118366$ |
| Reference (FeatFlow) | $5.57953523384$ | $0.010618948146$ | $0.11752016697$ |

## Grid convergence study (case 2D-1)

`benchmark.py` runs case 2D-1 for all five variants on a sequence of uniformly refined grids.
The advective term uses central differences (`Flux.UpwindWeight = 0.5`).

The sequence starts from two coarse Gmsh meshes. Each level is refined uniformly with `Grid.Refinement`
(each element is split into four), so the mesh size halves from level to level.

| Mesh | Mesh size channel / cylinder | Elements on level 0 |
| :--- | :---: | :--- |
| `meshes/cylinder_channel_coarse.geo` | $0.04$ / $0.004$ | 2416 triangles |
| `meshes/cylinder_channel_quad_coarse.geo` | $0.04$ / $0.004$ | 1236 quadrilaterals |

The FeatFlow reference values are only given for one (very fine) discretization. Therefore, the
errors are plotted over the total number of unknowns $N$ (velocity components and pressure):

* $N$ is the natural measure to compare discretizations with different dofs per element.
* The number of elements is listed in the tables as well.

The experimental order of convergence is computed with respect to $h \sim N^{-1/2}$:

$$
\text{EOC} = -2\,\frac{\log(e_l/e_{l-1})}{\log(N_l/N_{l-1})},
$$

where $e_l$ is the absolute error of $C_D$, $C_L$ or $\Delta p$ on level $l$.

Keep in mind that the refined grids keep the polygonal cylinder of level 0, which has 79 edges.
The geometric error of this polygon does not vanish under refinement. It eventually limits the
convergence of all indicators.

```bash
cd build-cmake/test/freeflow/navierstokes/flowaroundcylinder
make test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa \
     test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box \
     test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond \
     test_ff_navierstokes_dfg_benchmark_stationary_hybrid_pq1bubble_box
python3 benchmark.py                                       # all schemes, levels 0 to 3
python3 benchmark.py --schemes pq1bubble_box --levels 0 1  # a subset
python3 visualize.py
```

| Script / option | Purpose |
| :--- | :--- |
| `benchmark.py` | runs the simulations, writes `results/convergence_2d1.md` (errors and EOCs) |
| `--schemes`, `--levels` | subset of the variants and the refinement levels (default: all, 0 to 3) |
| `--memory-limit` | address space limit per run in GB (default 6), a run exceeding it fails |
| `--overwrite` | rerun runs that already have results |
| `-- <parameters>` | further parameters passed to the simulations |
| `visualize.py` | plots the indicators with the FeatFlow reference values and the bounds of the paper over $N$ to `results/convergence_2d1.png` |

In addition, each scheme runs once on the default mesh of the regression tests (`<scheme>_default`).
This mesh resolves the cylinder twice as fine as the coarse meshes (about 157 instead of 79 edges).
The plot shows these runs as unfilled markers, and the table lists them in a separate row.

Each run has its own folder `runs/2d1/<scheme>_r<level>/` (or `runs/2d1/<scheme>_default/`) in the build directory:

* the log file;
* the JSON file with the grid size, the number of unknowns and the indicators;
* the VTK output in `vtk/`.

A failed run is reported and skipped (missing data in the table and the plot).

## Running and evaluating the cases 2D-2 and 2D-3

`benchmark.py --case 2` and `--case 3` run one time-dependent case with the PQ1Bubble/Box variant, or with the
hybrid variant on triangles or quadrilaterals (`--scheme hybrid_pq1bubble_box` or `hybrid_pq1bubble_box_quad`),
each on its default mesh.
The runs take days on the default grid, so they are started in the background, detached from the terminal.

```bash
cd build-cmake/test/freeflow/navierstokes/flowaroundcylinder
python3 benchmark.py --case 2                        # 2D-2, default grid and time steps
python3 benchmark.py --case 3 --dt 0.00125 --tend 8  # 2D-3 with a larger time step
python3 benchmark.py --case 2 --scheme hybrid_pq1bubble_box  # 2D-2 with the hybrid variant
python3 visualize.py --case 2                        # evaluate (also while a run is going on)
python3 visualize.py --case 3
```

| Option | Purpose |
| :--- | :--- |
| `--scheme` | `pq1bubble_box` (default), `hybrid_pq1bubble_box` or `hybrid_pq1bubble_box_quad` |
| `--level` | `Grid.Refinement` of the default mesh (default 0) |
| `--dt` | fixed time step size (2D-2: the fine one after $3.5\,\text{s}$) |
| `--tend` | end time |
| `--time-scheme` | `TimeLoop.Scheme`: `DIRK2` (default), `DIRK3`, `ImplicitEuler` or `CrankNicolson` |
| `--foreground` | wait for the run instead of starting it in the background |
| `--overwrite` | replace the results of an existing run |
| `--tstart` (`visualize.py`) | 2D-2: the periodic regime is evaluated for $t \geq$ `tstart` (default $9\,\text{s}$, as the FeatFlow time series) |
| `--time-schemes` (`visualize.py`) | evaluate only the runs with these time stepping schemes (default `DIRK2`, `all` also includes older runs without the scheme in the name) |

Each run has its own folder `runs/2d<case>/<scheme>_r<level>[_dt<dt>]_<time scheme in lower case>/`, e.g.
`runs/2d2/pq1bubble_box_r0_dirk2/`:

* the log file and the process id (`<name>.pid`, stop the run with `kill $(cat <name>.pid)`);
* the CSV and JSON files with the time series of the indicators (the CSV file grows during the run);
* the VTK output in `vtk/`.

On the default grid (102341 unknowns, UMFPack, multithreaded assembly), one Newton iteration takes about
$2$ to $4\,\text{s}$ depending on the machine load. With DIRK2 (two stages with about three Newton iterations each),
one time step takes about $12$ to $22\,\text{s}$:

| Case | Time steps | Run time (DIRK2) |
| :--- | :--- | :--- |
| 2D-2 | 6535 | about one to one and a half days |
| 2D-3 | 12800 | about two to three days |

### Evaluation

`visualize.py` evaluates all runs of a case:

* **2D-2**: the last full period $[t_s, t_s + 1/f]$ after `--tstart`, where $t_s$ is a minimum of $C_L$
  (FeatFlow convention).
  * Max, min, mean $= (\max + \min)/2$ and amplitude $= \max - \min$ of $C_D$, $C_L$ and $\Delta p$.
  * The period $1/f$ and the Strouhal number $\mathrm{St} = D f/\bar{U}$.
  * $\Delta p(t_0 + 1/(2f))$, where $t_0$ is the time of the maximum of $C_L$ in the evaluated period, as
    defined in @cite SchaeferTurek1996. If $t_0 + 1/(2f)$ lies after the end of the run, the periodic value
    $\Delta p(t_0 - 1/(2f))$ is used.
  * Extrema of $C_L$ are refined with a parabola through three samples.
  * Minima of $C_L$ closer than $0.2\,\text{s}$ are merged, so noise near a minimum does not
    produce spurious periods.
* **2D-3**: the maxima of $C_D$ and $C_L$ and their times, and $\Delta p$ at $t = 8\,\text{s}$
  (the FeatFlow series end at the last half step $8 - \Delta t/2$, their value is extrapolated linearly).

The results go to three files in `results/`:

| File | Content |
| :--- | :--- |
| `case_2d_<case>.md` | all quantities of each run, compared with the bounds of @cite SchaeferTurek1996 and the finest FeatFlow $Q_2/P_1^\text{disc}$ run (level 6, $\Delta t = 1/1600\,\text{s}$) |
| `indicators_2d_<case>.png` | $C_{D,\max}$, $C_{L,\max}$ and $\Delta p$ of each run over its number of unknowns $N$, together with the FeatFlow results of all levels and time step sizes and the bounds of the paper |
| `case_2d_<case>.png` | the time series of each run (2D-2: $t \in [9, 10]\,\text{s}$, see below), plus `case_2d_2_olddata.png` with the old FeatFlow data |

All plots (also of case 2D-1) use the same conventions:

* The color shows the DuMux scheme.
* In the indicator plots, the marker shows the grid cell shape: triangles for the triangular grids, squares for the
  quadrilateral grids (also for the FeatFlow results). Several runs of the same scheme differ in the marker filling.
* The time series are plotted without markers, and several runs of the same scheme differ in the line style.

The time series are not given at the same points in time:

* FeatFlow stores the values at the half steps of the Crank-Nicolson method ($t_n + \Delta t/2$).
* DuMux stores them at the time steps $t_n$ (DIRK2; a Crank-Nicolson run also at the half steps).
* The time step sizes differ in case 2D-2 ($\Delta t = 0.001\,\text{s}$ vs. $1/1600\,\text{s}$).

Before plotting, `visualize.py` prints for each series how many of its time points coincide with the first FeatFlow
series, together with both time step sizes and the largest offset. For the plotted runs, the offset is at most
$\Delta t/2$, far below the period of about $0.33\,\text{s}$.

If the FeatFlow time series have been downloaded (see below), the time series plots also show FeatFlow results
on level 6:

* `case_2d_2.png`: the DuMux runs and FeatFlow $Q_2/P_1^\text{disc}$ (new data, $\Delta t = 1/1600\,\text{s}$),
  restricted to $t \in [9, 10]\,\text{s}$, the interval of the new data (the whole time range as long as no
  DuMux run has reached it).
* `case_2d_2_olddata.png`: the DuMux runs and the old data of both discretizations on the same level and time step size,
  $Q_2/P_1^\text{disc}$ and $\tilde{Q}_1/Q_0$ with $\Delta t = 1/800\,\text{s}$ (the old $Q_2/P_1^\text{disc}$ run
  with $1/1600\,\text{s}$ stops at $t = 0.025\,\text{s}$), on the same interval. The runs start from different
  states, so their phases differ.
* `case_2d_3.png`: the DuMux runs and FeatFlow $Q_2/P_1^\text{disc}$ ($\Delta t = 1/1600\,\text{s}$).
  FeatFlow provides no $\tilde{Q}_1/Q_0$ data for 2D-3.

### Videos

`video.py` renders a video of a run from its VTK output and the indicator time series (pyvista and matplotlib, MP4
by default, `--format gif` for a GIF):

```bash
python3 video.py --case 2 --run pq1bubble_box_r0_dirk2   # last three shedding periods
python3 video.py --case 3 --run pq1bubble_box_r0_dirk2   # whole run, t in [0, 8] s
```

* Panels: velocity magnitude, vorticity $\omega = \partial_x v - \partial_y u$ (vortex street), pressure, and
  $C_D$, $C_L$ over time with a marker at the current time.
* The fields are sampled on a uniform image of the channel (`--resolution`, pixels in $x$). The color ranges are
  fixed for the whole video: the 1st and 99th percentiles over up to 40 frames (symmetric for the vorticity).
* `--tstart`, `--tend` and `--stride` select the frames, `--fps` the frame rate. The video goes to
  `results/video_2d<case>_<run>.<format>` unless `--output` is given.

### Reference data

| File | Content |
| :--- | :--- |
| `reference/1996.json` | lower and upper bounds of @cite SchaeferTurek1996 per case |
| `reference/featflow.json` | FeatFlow results per case, written by `reference/fetch_featflow.py` |

`reference/featflow.json` holds the quantities of all FeatFlow runs of the cases 2D-2 and 2D-3. They are sorted by
the discretization, the time step size (coarse to fine) and the number of unknowns:

```
{"<case>": {"q2p1" | "q1p0": {"1/<steps per second>": [{"level", "numElements", "numDofs", "CDMax", ...}]}}}
```

* `q2p1` stands for $Q_2/P_1^\text{disc}$, `q1p0` for $\tilde{Q}_1/Q_0$ (rotated bilinear velocity, constant pressure).
* `numDofs` is the total number of scalar unknowns, as in the DuMux runs. The tables of the FeatFlow pages
  list a different number of degrees of freedom.

The FeatFlow time series (about $80\,\text{MB}$ with the archives) are not part of the repository.
`reference/fetch_featflow.py` handles them:

* It downloads the archives of the FeatFlow benchmark pages into `reference/featflow/<dataset>/` (ignored by git)
  and unpacks them there.
* It evaluates all runs with the same functions as the DuMux results (`visualize.py`).
* It writes the results to `reference/featflow.json`.

| Dataset | FeatFlow data | Used for |
| :--- | :--- | :--- |
| `2dcase2new` | 2D-2, $Q_2/P_1^\text{disc}$, levels 2 to 6, $\Delta t = 1/100, 1/400, 1/1600\,\text{s}$, $t \in [9, 10]\,\text{s}$ ("new reference results", higher-order integration of the forces) | 2D-2, `q2p1` |
| `2dcase2old` | 2D-2, $\tilde{Q}_1/Q_0$ and $Q_2/P_1^\text{disc}$, levels 2/3 to 6, $\Delta t = 1/200$ to $1/1600\,\text{s}$, $t \in [0, 10]\,\text{s}$ | 2D-2, `q1p0` (the old $Q_2/P_1^\text{disc}$ data is superseded) |
| `2dcase3` | 2D-3, $Q_2/P_1^\text{disc}$, levels 1 to 6, $\Delta t = 1/1600\,\text{s}$ | 2D-3, `q2p1` |

In the old data set, the forces are evaluated with a lower-order line integration. Compared with the new data set,
the drag is about $2\,\%$ too low on level 2 and $0.2\,\%$ on level 6 (pressure and period are identical).
Therefore, `indicators_2d_2.png` marks the $\tilde{Q}_1/Q_0$ results as old data.

Without DuMux runs, `visualize.py --case 2` and `--case 3` plot the references only.

```bash
python3 reference/fetch_featflow.py               # download (once), unpack and regenerate featflow.json
python3 reference/fetch_featflow.py --no-download # regenerate from the unpacked data
```

## Known limitations

* The time-dependent cases 2D-2 and 2D-3 are only implemented for the three PQ1Bubble/Box variants.
* The Crank-Nicolson method is not suitable for the pressure-dependent indicators (see above).
* The polygonal approximation of the cylinder limits the accuracy of all indicators on coarse grids.

## References
* FeatFlow benchmark suite, TU Dortmund:
  [2D-1 (Re = 20)](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html),
  [2D-2 (Re = 100, periodic)](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark2_re100.html),
  [2D-3 (Re = 100, fixed time interval)](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark3_re100.html).
