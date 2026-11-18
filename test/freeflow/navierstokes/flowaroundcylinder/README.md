# DFG Flow Around a Cylinder (Schäfer & Turek) {#benchmark-dfg-flow-around-cylinder}

This test implements the laminar flow around a cylinder benchmark of the DFG priority research
program "Flow Simulation with High-Performance Computers" @cite SchaeferTurek1996.

The benchmark defines a two-dimensional channel flow around a circular obstacle. Highly accurate
reference values are available for three quantities: the drag coefficient, the lift coefficient
and a pressure difference.

Currently, the stationary case 2D-1 ($\text{Re} = 20$) is implemented.

## Problem description

The incompressible Navier-Stokes equations

$$
\begin{aligned}
\rho\, \nabla\cdot(\mathbf{v}\otimes\mathbf{v}) - \mu\Delta\mathbf{v} + \nabla p &= \mathbf{0}, \\
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
| Inflow peak velocity | $U$ | $0.3$ | $\text{m/s}$ |
| Mean inflow velocity | $\bar{U} = \tfrac{2}{3}U$ | $0.2$ | $\text{m/s}$ |
| Reynolds number | $\text{Re} = \bar{U}D/\nu$ | $20$ | $-$ |
| Pressure evaluation point (cylinder front) | $A$ | $(0.15, 0.2)$ | $\text{m}$ |
| Pressure evaluation point (cylinder back) | $B$ | $(0.25, 0.2)$ | $\text{m}$ |

### Boundary conditions

* **Inflow** ($x = 0$): parabolic velocity profile $\mathbf{v}(0, y) = \left(4U y (H - y)/H^2,\ 0\right)^T$.
* **Walls** ($y = 0$, $y = H$) and **cylinder**: no-slip, $\mathbf{v} = \mathbf{0}$.
* **Outflow** ($x = 2.2$): do-nothing condition.
  * The advective momentum flux $\rho(\mathbf{v}\cdot\mathbf{n})\mathbf{v}$ leaves the domain.
  * The viscous and pressure contributions of the normal stress are zero.

## Quantities of interest (benchmark indicators)

The drag and lift forces act on the cylinder surface $S$. Let $\mathbf{n}$ be the unit normal on $S$
pointing into the cylinder and $v_\tau$ the tangential velocity. Then

$$
F_D = \int_S \left(\rho\nu\frac{\partial v_\tau}{\partial \mathbf{n}} n_y - p\, n_x\right) \text{d}S, \qquad
F_L = -\int_S \left(\rho\nu\frac{\partial v_\tau}{\partial \mathbf{n}} n_x + p\, n_y\right) \text{d}S .
$$

The benchmark indicators are

| Quantity | Definition | Reference value (FeatFlow) |
| :--- | :--- | :--- |
| Drag coefficient | $C_D = 2 F_D/(\rho\bar{U}^2 D)$ | $5.57953523384$ |
| Lift coefficient | $C_L = 2 F_L/(\rho\bar{U}^2 D)$ | $0.010618948146$ |
| Pressure difference | $\Delta p = p(A) - p(B)$ | $0.11752016697$ |

The reference values are taken from the FeatFlow benchmark page (see the references below).

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

The advective term is discretized with central differences (`Flux.UpwindWeight = 0.5`).

The unstructured grids are generated with Gmsh. Both use a mesh size of $h = 0.01$ in the channel and
$h = 0.002$ on the cylinder surface.

| Mesh | Elements |
| :--- | :--- |
| `meshes/cylinder_channel.geo` | 28948 triangles |
| `meshes/cylinder_channel_quad.geo` | 13826 quadrilaterals |

## Running the tests

```bash
cd build-cmake/test/freeflow/navierstokes/flow_around_cylinder
make test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box
./test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params_1.input
```

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

## Known limitations

* Only the stationary case 2D-1 is implemented. The time-dependent cases 2D-2 and 2D-3 ($\text{Re} = 100$)
  are not.
* The polygonal approximation of the cylinder limits the accuracy of all indicators on coarse grids.

## References
* FeatFlow benchmark suite, TU Dortmund:
  [2D-1 (Re = 20)](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html).
