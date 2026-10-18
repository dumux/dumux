# Shallow Water Poiseuille Flow {#benchmark-shallow-water-poiseuille}

## Steady channel flow with viscous friction at the side walls or at the bed

**Problem Description**

A straight channel of length $L$ and width $W$ has a bed that descends along the channel with the
slope $S$. A steady flow of uniform depth $h_0$ runs down the channel, so the free surface has the slope
of the bed. The water depth $h$ and the depth-averaged velocity components $u$ along and $v$ across the
channel satisfy the shallow water equations with the bed slope source term and a friction term. The
benchmark has two variants, which differ in the friction that balances gravity.

*Friction at the side walls.* The momentum balance contains the horizontal viscous flux
$\nabla \cdot (\nu_t h \nabla u)$ with a constant eddy viscosity $\nu_t$, and the velocity vanishes at the
side walls $y = \pm W/2$. For $h = h_0$, $v = 0$ and $u = u(y)$, the momentum balance along the channel
reduces to $0 = g h_0 S + \nu_t h_0 \, \mathrm{d}^2 u / \mathrm{d} y^2$, where $g$ is the gravitational
acceleration. Its solution is the parabolic profile

$$
u(y) = \frac{g S}{8 \nu_t} \left( W^2 - 4 y^2 \right).
$$

*Friction at the bed.* A laminar film flows down an inclined plane, with no slip at the bed and no shear
stress at the free surface. Its velocity profile over the depth is a half parabola with the mean velocity
$u$, and the bed shear stress is $\tau_b = 3 \mu u / h$ (@ref Dumux::FrictionLawViscousNoSlip), where $\mu$ is
the dynamic viscosity and $\rho$ the density of water. The side walls are free of friction, and the
depth-averaged momentum balance reduces to $0 = g h_0 S - \tau_b / \rho$, with the solution

$$
u = \frac{\rho g S h_0^2}{3 \mu}, \qquad v = 0,
$$

which is uniform in the channel.

**Parameters**

| Parameter                       | Symbol    | Side walls        | Bed                         | Unit              |
|---------------------------------|-----------|-------------------|-----------------------------|-------------------|
| Channel length                  | $L$       | $400$             | $1$                         | m                 |
| Channel width                   | $W$       | $100$             | $0.1$                       | m                 |
| Bed slope                       | $S$       | $5 \cdot 10^{-5}$ | $5 \cdot 10^{-6}$           | -                 |
| Water depth                     | $h_0$     | $10$              | $10^{-3}$                   | m                 |
| Eddy viscosity                  | $\nu_t$   | $0.1$             | -                           | m<sup>2</sup>/s   |
| Dynamic viscosity               | $\mu$     | -                 | $10^{-3}$                   | Pa s              |
| Density                         | $\rho$    | -                 | $1000$                      | kg/m<sup>3</sup>  |
| Gravitational acceleration      | $g$       | $9.81$            | $9.81$                      | m/s<sup>2</sup>   |

The maximum velocity of the flow with friction at the side walls is $u(0) = 6.131$ m/s. The mean velocity
of the film is $u = 1.635 \cdot 10^{-5}$ m/s.

**Setup**

The @ref ShallowWaterModels "shallow water model" with the viscous flux (`ShallowWater.EnableViscousFlux`)
is discretized with cell-centered finite volumes (@ref CCTpfaDiscretization). The flux over a face is the
flux of the exact solution of the Riemann problem with hydrostatic reconstruction of the bed @cite Audusse2004,
plus a viscous flux computed from the velocity difference between the adjacent cell centers. The bed
elevation is constant in each cell. Since the outside state at a boundary face has the bed of the inside
cell, a depth prescribed at a boundary face is shifted by the bed difference between face and cell center,
such that the free surface elevation at the face is matched.

*Friction at the side walls* (this folder, grid of $20 \times 20$ cells in `params.input`): the discharge
per unit width of the exact solution, $h_0 u(y)$, is prescribed at the inflow and the depth $h_0$ at the
outflow. At the side walls, the viscous flux is computed from the velocity difference between the cell
center and the wall (`Problem.AlphaWall = 1`, no slip). The simulation starts from the exact solution and
uses the implicit Euler method up to $t = 2880$ s, when the flow is steady: on the grid of `params.input`, the
errors at $t = 2880$ s and $t = 10^6$ s agree in all printed digits. The tests `test_shallowwater_poiseuilleflow`
(and `_parallel`, `_unstructured` on a triangular grid, `_unstructured_parallel`) compare the solution with
a reference solution.

*Friction at the bed* (folder `vertical`, grid of $20 \times 2$ cells): the depth $h_0$ is prescribed at
inflow and outflow, and the side walls are free of friction. The bed shear stress enters the momentum
balance as a source term. The stationary problem is solved directly, starting from the exact solution. The
test `test_shallowwater_poiseuilleflow_vertical` compares the solution with a reference solution.

**Results**

The figure and the tables are produced by running
```bash
python3 plot_results.py ./test_shallowwater_poiseuilleflow vertical/test_shallowwater_poiseuilleflow_vertical
```
in the build directory `test/freeflow/shallowwater/poiseuilleflow` (requires NumPy and Matplotlib).
The errors are the discrete $L^2$ norms $\| e \|_2 = \left( \sum_i |V_i| \, e_i^2 / \sum_i |V_i| \right)^{1/2}$,
where $e_i$ is the difference to the exact velocity at the center of the cell $V_i$ with area $|V_i|$.
The flow with friction at the side walls is computed up to $t = 10^6$ s with time steps of up to
$10^5$ s.

![Error of the velocity profile in the middle of the channel (left), errors of the flow with friction at the side walls (center) and of the film flow (right) under grid refinement](shallowwater_poiseuille.png){html: width=100%}

*Friction at the side walls.* The grid is refined along the channel with 320 cells across, and across
the channel with 320 cells along.

| grid (along × across) | $L^2$ error $u$ (m/s) | rate |
|--:|--:|--:|
| 10 × 320 | 7.954e-04 | - |
| 20 × 320 | 3.902e-04 | 1.03 |
| 40 × 320 | 2.045e-04 | 0.93 |
| 80 × 320 | 1.130e-04 | 0.86 |

| grid (along × across) | $L^2$ error $u$ (m/s) | rate |
|--:|--:|--:|
| 320 × 10 | 6.320e-03 | - |
| 320 × 20 | 2.741e-03 | 1.21 |
| 320 × 40 | 7.782e-04 | 1.82 |
| 320 × 80 | 2.106e-04 | 1.89 |

Along the channel the error converges with first order, across the channel with second order. The rates
along the channel decrease as the error approaches that of the 320 cells across; with 160 cells across, the
convergence test below obtains the rates 1.01 and 0.88. Across the channel, the largest errors are in the
cells next to the walls.

*Friction at the bed.* The grid is refined along the channel.

| grid (along × across) | $L^2$ error $u$ (m/s) | rate |
|--:|--:|--:|
| 20 × 2 | 1.947e-09 | - |
| 40 × 2 | 1.001e-09 | 0.96 |
| 80 × 2 | 5.063e-10 | 0.98 |
| 160 × 2 | 2.545e-10 | 0.99 |

The error converges with first order. On 20 cells, the relative error of the velocity is
$1.19 \cdot 10^{-4}$, close to $S \Delta x / (2 h_0) = 1.25 \cdot 10^{-4}$, where $\Delta x$ is the cell
size along the channel and $S \Delta x$ the step of the bed at a face. The script repeats the 20-cell run with
the slope $10 S$ and with the depth $2 h_0$; the relative errors are $1.19 \cdot 10^{-3}$ and
$6.04 \cdot 10^{-5}$, again 0.95 to 0.97 times $S \Delta x / (2 h_0)$.

The test `test_shallowwater_poiseuilleflow_convergence` refines the grid along the channel
($10, 20, 40$ cells, 160 cells across) and across the channel ($20, 40, 80$ cells, 160 cells along) and
checks that the mean convergence rates are at least 0.85 and 1.6. The test
`test_shallowwater_poiseuilleflow_vertical_convergence` runs the film flow on $20, 40, 80$ cells and checks
that the mean rate is at least 0.9.
