# Oscillation in a Paraboloid Bowl {#benchmark-thacker-bowl}

## Wetting and drying with the exact solution of Thacker

**Problem Description**

Water oscillates without friction in a bowl with the bed elevation $z(r) = D_0 r^2 / L^2$, where
$r = \sqrt{x^2 + y^2}$ is the distance from the center of the bowl. When the water is at rest, its free
surface is at the elevation $D_0$ and meets the bed at the radius $L$. Initially, the water is at rest with
a free surface that is a paraboloid raised by $\eta_0$ at the center. The free surface oscillates, and the
shoreline, where the water depth $h$ vanishes, moves up and down the bed of the bowl. The shallow water
equations hold in the wet part of the domain $[-2\,\mathrm{m}, 2\,\mathrm{m}]^2$.

**Exact Solution**

Thacker @cite Thacker1981 (Section 4) gives the periodic, radially symmetric solution (see also the
compilation of shallow water solutions SWASHES @cite Delestre2013). With the frequency
$\omega = \sqrt{8 g D_0}/L$, the gravitational acceleration $g$, the amplitude parameter
$A = \left( (D_0 + \eta_0)^2 - D_0^2 \right) / \left( (D_0 + \eta_0)^2 + D_0^2 \right)$ and
$c(t) = 1 - A \cos \omega t$, the water depth and the velocity components in $x$- and $y$-direction are

$$
h(r, t) = \max\left( D_0 \left( \frac{\sqrt{1 - A^2}}{c(t)} - \frac{r^2}{L^2} \frac{1 - A^2}{c(t)^2} \right), 0 \right),
\qquad
(u, v) = \frac{\omega A \sin \omega t}{2 c(t)} (x, y)
$$

in the wet region $h > 0$, and $u = v = 0$ outside of it. The free surface remains a paraboloid. The period
of the oscillation is $T = 2 \pi / \omega$. For the parameters below, $T = 6.5518$ s and $A = 0.15294$.

**Parameters**

| Parameter                                | Symbol           | Value     | Unit             |
|------------------------------------------|------------------|-----------|------------------|
| Depth of the bowl below the water at rest | $D_0$           | $0.03$    | m                |
| Radius of the shoreline at rest          | $L$              | $1.6$     | m                |
| Initial elevation of the free surface at the center | $\eta_0$ | $0.005$ | m              |
| Gravitational acceleration               | $g$              | $9.81$    | m/s<sup>2</sup>  |
| Final time                               | $t_\mathrm{end}$ | $2 T$     | s                |

**Setup**

The @ref ShallowWaterModels "shallow water model" is discretized with cell-centered finite volumes
(@ref CCTpfaDiscretization) and the implicit Euler method. The flux over a face is the flux of the exact
solution of the Riemann problem with hydrostatic reconstruction of the bed @cite Audusse2004. The bed
elevation is constant in each cell and is evaluated at the cell center. The boundaries are walls.
Where the bed is dry at $t = 0$, the initial depth is a film of $10^{-5}$ m. The mass flux over a face is
multiplied by a factor that decreases from one to zero as the upwind depth falls below
$10^{-5}$ m (`FluxLimiterLET`).

The grid of `params.input` has $40 \times 40$ cells and the time step size is $\Delta t = 0.1$ s. The
test `test_shallowwater_bowl` (and `_parallel` on two processes) compares the solution at $t = 13$ s with
a reference solution. The test prints the errors with respect to the exact solution after each time step and
writes the exact solution to the output (`exactWaterDepth`, `exactVelocityX`, `exactVelocityY`). For the
grid study, the grid has $N \times N$ cells and the time step size is $\Delta t = 4\,\mathrm{s}/N$.

**Results**

The figure and the table are produced by running
```bash
python3 plot_results.py ./test_shallowwater_bowl
```
in the build directory `test/freeflow/shallowwater/bowl` (requires NumPy and Matplotlib). The script
computes the exact solution independently of DuMux.

![Free surface (left) and discharge (center) along the row of cells next to y = 0, and L1 error of the water depth over time (right)](shallowwater_bowl.png){html: width=100%}

The table lists the discrete $L^1$ errors of the water depth,
$\| e \|_1 = \sum_i |V_i| \, |e_i| / \sum_i |V_i|$, where $e_i$ is the difference to the exact solution
at the center of the cell $V_i$ with area $|V_i|$, at the time step closest to $T$ and at $t_\mathrm{end}$.

| grid | $\Delta t$ (s) | $L^1$ error $h$, $t = T$ (m) | rate | $L^1$ error $h$, $t = 2 T$ (m) | rate |
|--:|--:|--:|--:|--:|--:|
| 40 × 40 | 0.1000 | 8.736e-04 | - | 1.149e-03 | - |
| 80 × 80 | 0.0500 | 5.674e-04 | 0.62 | 9.006e-04 | 0.35 |
| 160 × 160 | 0.0250 | 3.218e-04 | 0.82 | 5.591e-04 | 0.69 |
| 320 × 320 | 0.0125 | 1.701e-04 | 0.92 | 2.951e-04 | 0.92 |

The rates increase towards one under refinement. The computed oscillation is damped: at $t = 11$ s, close to
the largest velocity of the second period, the largest discharge $h u$ in the row of cells next to $y = 0$ is
0.041, 0.209, 0.487 and 0.728 times the exact one on the four grids. An eight times smaller time step
on the coarsest grid raises the ratio from 0.041 to 0.059, a refinement of the grid with half the time step
to 0.209: the damping comes mainly from the spatial discretization.
