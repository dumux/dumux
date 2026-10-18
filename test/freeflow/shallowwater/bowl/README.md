# Oscillation in a Paraboloid Bowl {#benchmark-thacker-bowl}

## Wetting and drying with the exact solution of Thacker

**Problem Description**

Water oscillates without friction in a bowl with the bed elevation $z(r) = D_0 r^2 / L^2$, where
$r = \sqrt{x^2 + y^2}$ is the distance from the center of the bowl. When the water is at rest, its free
surface is at the elevation $D_0$ and meets the bed at the radius $L$. Initially, the water is at rest with
a free surface that is a paraboloid raised by $\eta_0$ at the center. The free surface oscillates, and the
shoreline, where the water depth $h$ vanishes, moves up and down the bed of the bowl. The shallow water
equations hold in the wet part of the domain $[-2\,\mathrm{m}, 2\,\mathrm{m}]^2$.

![Setup of the bowl, cross section through its center with exaggerated vertical scale](shallowwater_bowl_setup.svg){html: width=70%}

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

The figure and the tables are produced by running
```bash
python3 plot_results.py ./test_shallowwater_bowl
```
in the build directory `test/freeflow/shallowwater/bowl` (requires NumPy and Matplotlib). The script
computes the exact solution independently of DuMux.

![Free surface and discharge along the row of cells next to y = 0, and relative L1 error of the water depth over time](shallowwater_bowl.png){html: width=100%}

The table lists the relative discrete $L^1$ errors of the water depth, $\| e \|_1 / \| h \|_1$ with
$\| e \|_1 = \sum_i |V_i| \, |e_i|$, where $e_i$ is the difference to the exact solution at the center of
the cell $V_i$ with area $|V_i|$, and $\| h \|_1$ the norm of the exact depth evaluated in the same way, at the
time step closest to $T$ and at $t_\mathrm{end}$. The $L^1$ norm is used because the solution is not smooth at
the shoreline, where the depth has a kink and the velocity jumps to zero.

| grid | $\Delta t$ (s) | rel. $L^1$ error $h$, $t = T$ | rate | rel. $L^1$ error $h$, $t = 2 T$ | rate |
|--:|--:|--:|--:|--:|--:|
| 40 × 40 | 0.1000 | 1.159e-01 | - | 1.524e-01 | - |
| 80 × 80 | 0.0500 | 7.525e-02 | 0.62 | 1.195e-01 | 0.35 |
| 160 × 160 | 0.0250 | 4.269e-02 | 0.82 | 7.415e-02 | 0.69 |
| 320 × 320 | 0.0125 | 2.256e-02 | 0.92 | 3.914e-02 | 0.92 |

On the finer grids the error converges with approximately first order; on the coarser grids, and more so
at $t = 2 T$ than at $t = T$, it converges more slowly. The reason is the damping of the oscillation. The next
table lists the largest discharge $h u$ in the row of cells next to $y = 0$ relative to the exact one at
$t = 11$ s, close to the largest velocity of the second period:

| grid | $\Delta t$ (s) | peak discharge / exact peak discharge | logarithm |
|--:|--:|--:|--:|
| 40 × 40 | 0.1000 | 0.041 | -3.21 |
| 80 × 80 | 0.0500 | 0.209 | -1.57 |
| 160 × 160 | 0.0250 | 0.487 | -0.72 |
| 320 × 320 | 0.0125 | 0.728 | -0.32 |
| 40 × 40 | 0.0125 | 0.059 | -2.83 |

The logarithm of the ratio approximately halves with each refinement: the oscillation decays exponentially,
with a decay rate approximately proportional to the cell size. On the coarse grids the oscillation has
almost decayed, and the error approaches the amplitude of the oscillation instead of decreasing in
proportion to the cell size. An eight times smaller time step on the coarsest grid raises the ratio from 0.041
to 0.059, a refinement of the grid with half the time step to 0.209: the damping comes mainly from the spatial
discretization.

The wet/dry front is not the cause of the slow convergence. The script splits the relative $L^1$ error
of the depth into the contributions of the cells within 0.2 m of the exact shoreline, of the wet interior
and of the dry bed:

| $t$ (s) | grid | near the shoreline | rate | wet interior | rate | dry bed | rate |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 7 | 40 × 40 | 3.875e-02 | - | 5.628e-02 | - | 5.467e-04 | - |
| 7 | 80 × 80 | 2.087e-02 | 0.89 | 3.649e-02 | 0.63 | 5.653e-04 | -0.05 |
| 7 | 160 × 160 | 8.716e-03 | 1.26 | 2.135e-02 | 0.77 | 5.604e-04 | 0.01 |
| 7 | 320 × 320 | 3.612e-03 | 1.27 | 1.113e-02 | 0.94 | 4.721e-04 | 0.25 |
| 13 | 40 × 40 | 6.419e-02 | - | 8.809e-02 | - | 5.144e-04 | - |
| 13 | 80 × 80 | 5.072e-02 | 0.34 | 7.022e-02 | 0.33 | 5.035e-04 | 0.03 |
| 13 | 160 × 160 | 3.124e-02 | 0.70 | 4.506e-02 | 0.64 | 4.003e-04 | 0.33 |
| 13 | 320 × 320 | 1.614e-02 | 0.95 | 2.477e-02 | 0.86 | 2.502e-04 | 0.68 |

The contribution of the wet interior is the larger one on all grids and converges like the total error. The
contribution of the dry bed is the water film of initially $10^{-5}$ m and stays small.
