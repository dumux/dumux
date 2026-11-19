# DFG Flow Around a Cylinder (Schäfer & Turek) {#benchmark-dfg-flow-around-cylinder}

This test implements the laminar flow around a cylinder benchmark of the DFG priority research
program "Flow Simulation with High-Performance Computers" @cite SchaeferTurek1996.

The benchmark defines a two-dimensional channel flow around a circular obstacle. Highly accurate
reference values are available for three quantities: the drag coefficient, the lift coefficient
and a pressure difference.

All three 2D cases are implemented:

* the stationary case 2D-1 ($\text{Re} = 20$) for all discretization variants;
* the time-dependent cases 2D-2 and 2D-3 ($\text{Re} = 100$) for the PQ1Bubble/Box variant.

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

For case 2D-2, the frequency $f$ of the vortex shedding defines the Strouhal number $\text{St} = D f/\bar{U}$.

The FeatFlow reference values are listed below (see the references at the end):

* 2D-2 and 2D-3: $Q_2/P_1^\text{disc}$ elements on level 6, Crank-Nicolson method with $\Delta t = 1/1600\,\text{s}$;
* 2D-2: the Strouhal number is taken from the FeatFlow benchmark page (level 6, $\Delta t = 1/200\,\text{s}$).

| Case | Quantity | Reference value |
| :--- | :--- | :--- |
| 2D-1 | $C_D$ | $5.57953523384$ |
| 2D-1 | $C_L$ | $0.010618948146$ |
| 2D-1 | $\Delta p$ | $0.11752016697$ |
| 2D-2 | $C_{D,\max}$ | $3.22739$ |
| 2D-2 | $C_{L,\max}$ | $0.98657$ |
| 2D-2 | $\text{St}$ | $0.30188$ |
| 2D-3 | $C_{D,\max}$ at $t$ | $2.9437637214$ at $t = 3.9365625\,\text{s}$ |
| 2D-3 | $C_{L,\max}$ at $t$ | $0.47748781595$ at $t = 5.6928125\,\text{s}$ |
| 2D-3 | $\Delta p(t = 8\,\text{s})$ | $-0.11154138872$ |

### Evaluation in the code

* **Drag and lift**: the traction is integrated over the boundary faces on the cylinder.
  On the no-slip cylinder wall, this is equivalent to the formula above.
* **Pressure difference**: the pressure is interpolated at $A$ and $B$.
  * The cylinder surface is approximated by a polygon, so $A$ and $B$ may lie outside the grid.
  * In that case, they are moved upstream ($A$) or downstream ($B$) in steps of $1\,\text{mm}$ until
    they hit an element.
* The evaluation is not parallelized. A parallel run with inertia terms stops with an exception.

## Discretization

The momentum balance and the mass balance are discretized separately and coupled by the
`FreeFlowCouplingManager`. Three combinations of discretization schemes are tested:

| Momentum | Mass | Grid |
| :--- | :--- | :--- |
| @ref DiamondDiscretization | @ref CCTpfaDiscretization | triangles |
| @ref PQ1BubbleDiscretization | @ref BoxDiscretization | triangles |
| @ref PQ1BubbleDiscretization | @ref DiamondDiscretization | quadrilaterals |

The PQ1Bubble/Box variant uses the new CVFE problem and variables interfaces (`problem_newinterface.hh`,
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

The time-dependent cases use the Crank-Nicolson method ($\theta$-scheme with $\theta = 0.5$), as the FeatFlow
reference computations.

* It is second-order accurate. First-order schemes such as the implicit Euler method damp the oscillations
  @cite SchaeferTurek1996.
* It is run with the multi-stage assembler and time stepper (`Experimental::MultiStageMultiDomainAssembler`,
  `Experimental::MultiStageTimeStepper`). These are only available for the CVFE schemes with the new problem
  interface, i.e. the PQ1Bubble/Box variant.
* The time-dependent inflow (Dirichlet constraints) is updated at each stage time.
* Newton's method uses the shift criterion (`Newton.MaxRelativeShift`), since the multi-stage assembler
  cannot assemble the residual alone.

Fixed time step sizes are used, following the benchmark description:

| Case | Time steps | Parameters |
| :--- | :--- | :--- |
| 2D-2 | coarse step $0.1\,\text{s}$ until $3.5\,\text{s}$ to develop the periodic flow, then $0.001\,\text{s}$ until $30\,\text{s}$ | `TimeLoop.Dt1`, `T1`, `Dt2`, `T2` |
| 2D-3 | $1/1600\,\text{s}$ on $[0, 8]\,\text{s}$ | `TimeLoop.Dt`, `TEnd` |

If Newton's method fails, the step is repeated with a reduced step size. Afterwards, the fixed step size is used again.

### Evaluation at the half step

With $\theta = 0.5$, the weighting also applies to the pressure, which is the Lagrange multiplier of the
incompressibility constraint, and to the continuity equation. Therefore:

* Only the average $\tfrac{1}{2}(p^n + p^{n+1})$ is determined. It approximates the pressure at $t^{n+1/2}$.
* The pressure at the time levels $t^n$ shows an undamped odd-even oscillation, excited by the inconsistent
  initial pressure $p^0 = 0$. The VTK output (time levels) contains this oscillation.
* The benchmark indicators are evaluated from the averaged solution $\tfrac{1}{2}(\mathbf{x}^n + \mathbf{x}^{n+1})$
  at $t^{n+1/2}$. They are free of the oscillation, and the FeatFlow reference time series are given at
  the half steps as well.

The indicators of every time step are written to the JSON file `<Problem.Name>.json` (arrays `time`, `cDrag`,
`cLift`, `pDiff`). The file is written at the end of the simulation and also if it fails (`finished`, `error`).

## Running the tests

```bash
cd build-cmake/test/freeflow/navierstokes/flow_around_cylinder
make test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box
./test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params_1.input  # case 2D-1
./test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params_2.input  # case 2D-2 (days)
./test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params_3.input  # case 2D-3 (days)
```

The ctest regression tests only cover the stationary case 2D-1.

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
| PQ1Bubble / Diamond (quadrilaterals) | $5.61234$ | $0.00387379$ | $0.118366$ |
| Reference (FeatFlow) | $5.57953523384$ | $0.010618948146$ | $0.11752016697$ |

## Grid convergence study (case 2D-1)

`benchmark.py` runs case 2D-1 for all three variants on a sequence of uniformly refined grids.
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
cd build-cmake/test/freeflow/navierstokes/flow_around_cylinder
make test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa \
     test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box \
     test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond
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

Each run has its own folder `runs/2d1/<scheme>_r<level>/` in the build directory:

* the log file;
* the JSON file with the grid size, the number of unknowns and the indicators;
* the VTK output in `vtk/`.

A failed run is reported and skipped (missing data in the table and the plot).

## Known limitations

* The time-dependent cases 2D-2 and 2D-3 are only implemented for the PQ1Bubble/Box variant.
* The Crank-Nicolson method is not L-stable. The impulsive start of case 2D-2 from rest may cause
  oscillations in the beginning, which are not damped.
* The polygonal approximation of the cylinder limits the accuracy of all indicators on coarse grids.

## References
* FeatFlow benchmark suite, TU Dortmund:
  [2D-1 (Re = 20)](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html),
  [2D-2 (Re = 100, periodic)](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark2_re100.html),
  [2D-3 (Re = 100, fixed time interval)](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark3_re100.html).
