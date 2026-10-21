# Wet Dam Break {#benchmark-dam-break}

## Riemann problem of the shallow water equations on a horizontal bed

**Problem Description**

A horizontal channel without friction, of length $L$ and width $W$, is divided by a gate at $x = x_0$.
At $t = 0$ the water is at rest, with the depth $h_L$ upstream and $h_R < h_L$ downstream of the gate,
and the gate is removed at once. All boundaries are walls.

![Setup of the dam break](shallowwater_dambreak_setup.svg){html: width=70%}

The water depth $h$ and the depth-averaged velocity $u$ along the channel satisfy the one-dimensional
shallow water equations

$$
\frac{\partial h}{\partial t} + \frac{\partial (h u)}{\partial x} = 0, \qquad
\frac{\partial (h u)}{\partial t} + \frac{\partial}{\partial x} \left( h u^2 + \frac{1}{2} g h^2 \right) = 0,
$$

where $g$ is the gravitational acceleration. The flow does not vary across the channel and the
velocity component $v$ across the channel vanishes.

**Exact Solution**

Until the waves reach the end walls, the solution is the self-similar solution of the Riemann problem
with the states $(h_L, 0)$ and $(h_R, 0)$ given by Stoker @cite Stoker1957 (see also the compilation
of shallow water solutions SWASHES @cite Delestre2013). A rarefaction wave travels upstream and a bore
downstream. Between them, the depth $h_m$ and the velocity $u_m$ are constant. With the celerities
$c_L = \sqrt{g h_L}$ and $c_m = \sqrt{g h_m}$, the depth $h_m$ is the root of

$$
2 \left( c_L - c_m \right) = \left( h_m - h_R \right) \sqrt{\frac{g}{2} \frac{h_m + h_R}{h_m h_R}},
$$

where the left side is the velocity behind the rarefaction wave and the right side the velocity behind
the bore. Then $u_m = 2 (c_L - c_m)$ and the bore moves with the speed $s = h_m u_m / (h_m - h_R)$.
With $\xi = (x - x_0)/t$, the solution is

$$
(h, u) = \begin{cases}
(h_L, 0), & \xi \leq -c_L, \\
\left( (2 c_L - \xi)^2 / (9 g), \; 2 (c_L + \xi) / 3 \right), & -c_L < \xi \leq u_m - c_m, \\
(h_m, u_m), & u_m - c_m < \xi \leq s, \\
(h_R, 0), & \xi > s.
\end{cases}
$$

For the parameters below, $h_m = 2.2070$ m, $u_m = 3.2223$ m/s and $s = 5.8921$ m/s. At the final time
$t_\mathrm{end}$, the head of the rarefaction wave is at $x = 3.74$ m and the bore at $x = 15.89$ m.

**Parameters**

| Parameter                       | Symbol             | Value  | Unit             |
|---------------------------------|--------------------|--------|------------------|
| Channel length                  | $L$                | $20$   | m                |
| Channel width                   | $W$                | $5$    | m                |
| Gate position                   | $x_0$              | $10$   | m                |
| Initial depth upstream          | $h_L$              | $4$    | m                |
| Initial depth downstream        | $h_R$              | $1$    | m                |
| Gravitational acceleration      | $g$                | $9.81$ | m/s<sup>2</sup>  |
| Final time                      | $t_\mathrm{end}$   | $1$    | s                |

**Setup**

The @ref ShallowWaterModels "shallow water model" is discretized with cell-centered finite volumes
(@ref CCTpfaDiscretization) and the implicit Euler method. The flux over a face is the flux of the exact
solution of the Riemann problem between the states of the two adjacent cells. At the walls, the outside
state is the inside state with the normal velocity reversed. The scheme is first order in space and time.
The test computes the exact solution with the same Riemann solver and writes it to the output
(`exactWaterDepth`, `exactVelocityX`).

The grid of `params.input` has 100 cells along and 25 cells across the channel, and the time step size is
$\Delta t = 0.005$ s. The test `test_shallowwater_dambreak_sequential` (and `_parallel` on two processes)
compares the solution at $t_\mathrm{end}$ with a reference solution. For the grid study, the channel is
resolved with $N$ cells along and one cell across the channel, and the time step size is
$\Delta t = 0.5\,\mathrm{s}/N$, such that the Courant number $c_L \Delta t / \Delta x = 0.16$, with the
cell size $\Delta x = L/N$, is the same on all grids.

*Second-order flux.* With the type tag `ShallowWaterSecondOrder` (@ref Dumux::ShallowWaterSecondOrderFlux),
the states of both cells are reconstructed linearly at each face before the Riemann problem is solved. Water
depth, free surface elevation and velocity get least-squares gradients over the face neighbors, limited at each
face with a smooth van Albada average (@ref Dumux::ShallowWater::CellReconstruction). The test
`test_shallowwater_dambreak_secondorder_convergence` runs the grids with 100, 200 and 400 cells with this flux,
checks the same minimum convergence rate, and checks that the water volume is conserved to $10^{-12}$.

**Results**

The figures and the tables are produced by running
```bash
python3 plot_results.py ./test_shallowwater_dambreak --second-order ./test_shallowwater_dambreak_secondorder
```
in the build directory `test/freeflow/shallowwater/dambreak` (requires NumPy, SciPy and Matplotlib).
The script computes the exact solution independently of DuMux. It agrees with the exact solution written by
the test to within $10^{-13}$ m. For comparison, the script also solves the dam break with a textbook
first-order Godunov scheme, the exact Riemann solver with the explicit Euler method, on the same grids and
with the same time step sizes.

![Water depth and velocity at t = 1 s, and relative L1 errors under grid refinement](shallowwater_dambreak.png){html: width=100%}

The table lists the relative discrete $L^1$ errors at $t_\mathrm{end}$,
$\| e \|_1 / \| h \|_1$ for the depth and $\| e \|_1 / \| u \|_1$ for the velocity, with
$\| e \|_1 = \sum_i |V_i| \, |e_i|$, where $e_i$ is the difference to the exact solution at the center of
the cell $V_i$ with area $|V_i|$, and the norms of the exact solution evaluated in the same way. The $L^1$
norm is used because the solution is discontinuous at the bore: there, the error is of the order of the jump
over a width proportional to $\Delta x$, which gives $L^1$ errors proportional to $\Delta x$ but $L^2$ errors
proportional to $\Delta x^{1/2}$.

| $N$ | $\Delta t$ (s) | rel. $L^1$ error $h$ | rate | rel. $L^1$ error $u$ | rate | rel. $L^1$ error $h$, explicit Godunov | rate |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 100 | 5.000e-03 | 2.841e-02 | - | 9.613e-02 | - | 2.330e-02 | - |
| 200 | 2.500e-03 | 1.633e-02 | 0.80 | 5.346e-02 | 0.85 | 1.325e-02 | 0.81 |
| 400 | 1.250e-03 | 9.405e-03 | 0.80 | 3.047e-02 | 0.81 | 7.631e-03 | 0.80 |
| 800 | 6.250e-04 | 5.392e-03 | 0.80 | 1.727e-02 | 0.82 | 4.381e-03 | 0.80 |
| 1600 | 3.125e-04 | 3.063e-03 | 0.82 | 9.826e-03 | 0.81 | 2.475e-03 | 0.82 |

Both errors converge with approximately first order. The explicit Godunov scheme converges with the same
rates, and its errors are smaller by a factor of 1.22 to 1.24. The rates are therefore a property of
first-order schemes for this solution and not of the implicit time discretization. The script splits the
relative $L^1$ error of the depth into the contributions of the cells within 1 m of each wave. Near the
bore, the contribution converges with approximately first order, from $9.93 \times 10^{-3}$ on 100 cells to
$6.45 \times 10^{-4}$ on 1600 cells. Near the head of the rarefaction wave, where the exact depth has a kink,
it converges more slowly, from $6.04 \times 10^{-3}$ to $1.08 \times 10^{-3}$, and makes up 35 % of the error
on 1600 cells.

The test `test_shallowwater_dambreak_convergence` runs the grids with 100, 200 and 400 cells and checks
that the mean convergence rates of the relative $L^1$ errors of depth and velocity are at least 0.7.

*Second-order flux.* The script runs the same grids and time step sizes with the second-order flux.

![Water depth and velocity at t = 1 s on 400 cells with the first- and the second-order flux, and relative L1 errors of both under grid refinement](shallowwater_dambreak_secondorder.png){html: width=100%}

| $N$ | rel. $L^1$ error $h$ | rate | rel. $L^1$ error $u$ | rate | error $h$ / first order | $\max h - h_L$ (m) | $h_R - \min h$ (m) |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 100 | 9.649e-03 | - | 3.527e-02 | - | 0.34 | 1.3e-03 | 1.6e-04 |
| 200 | 4.733e-03 | 1.03 | 1.614e-02 | 1.13 | 0.29 | 1.2e-03 | 1.4e-04 |
| 400 | 2.526e-03 | 0.91 | 8.474e-03 | 0.93 | 0.27 | 9.0e-04 | 1.3e-04 |
| 800 | 1.364e-03 | 0.89 | 4.477e-03 | 0.92 | 0.25 | 4.9e-04 | 1.4e-04 |
| 1600 | 7.523e-04 | 0.86 | 2.569e-03 | 0.80 | 0.25 | 1.8e-04 | 1.7e-04 |

The errors of the depth are 25 % to 34 % of those of the first-order flux, and both errors converge with
approximately first order, as the discontinuity of the bore permits. In the profiles on 400 cells, the head of
the rarefaction wave is resolved more sharply. The depth leaves the range $[h_R, h_L]$ of the exact solution by
at most $1.3 \cdot 10^{-3}$ m. These over- and undershoots are smaller than the smoothing scale of the limiter,
$10^{-3}$ times the depth of the cell (`ShallowWater.LimiterRelativeEpsilon`), below which increments are not
limited.
