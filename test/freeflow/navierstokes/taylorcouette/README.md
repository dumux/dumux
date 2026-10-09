# Taylor-Couette {#benchmark-taylor-couette}

## Steady viscous flow between two rotating cylinders

**Problem Description**

We simulate the steady, laminar flow of an incompressible fluid confined in the annular gap between two concentric, independently rotating cylinders — the classical Taylor-Couette configuration @cite taylor1923. Because a closed-form solution exists, this test case allows a direct comparison of the numerical velocity and pressure fields against the analytical solution.

Assumptions:

- steady state
- laminar flow
- incompressible fluid
- no gravity

The inner cylinder of radius $R_1$ rotates at angular velocity $\Omega_1$, and the outer cylinder of radius $R_2$ rotates at angular velocity $\Omega_2$. In the reference configuration, the outer cylinder is fixed ($\Omega_2 = 0$) and the inner cylinder rotates at $\Omega_1 = 100$ rad/s, corresponding to a Reynolds number $Re = |\Omega_1 R_1| (R_2 - R_1) / \nu = 100$.

**Analytical Reference Solution**

Taylor @cite taylor1923 showed that the tangential velocity $u_{\theta}$ within the annular gap is described by:

$$
u_{\theta}(r) = A r + \frac{B}{r}, \qquad
A = \frac{\Omega_2 R_2^2 - \Omega_1 R_1^2}{R_2^2 - R_1^2},
\qquad
B = \frac{(\Omega_1 - \Omega_2) R_1^2 R_2^2}{R_2^2 - R_1^2}
$$

The corresponding pressure field follows from the radial momentum balance $\partial p/\partial r = \rho\, u_\theta^2/r$:

$$
p(r) = \frac{A^2 r^2}{2} + 2AB \ln(r) - \frac{B^2}{2r^2} + C
$$

where $C$ is fixed by a reference pressure condition.

**Parameters**

| Parameter                        | Symbol     | Value  | Unit  |
|-----------------------------------|------------|--------|-------|
| Inner cylinder radius              | $R_1$      | $1$    | m     |
| Outer cylinder radius              | $R_2$      | $2$    | m     |
| Inner cylinder angular velocity    | $\Omega_1$ | $100$  | rad/s |
| Outer cylinder angular velocity    | $\Omega_2$ | $0$    | rad/s |
| Fluid density                      | $\rho$     | $1$    | kg/m³ |
| Fluid kinematic viscosity          | $\nu$      | $1$    | m²/s  |

**Setup**

The implementation uses the DuMux free-flow Navier-Stokes model. The annular domain is discretized using DuMux's [`CakeGridManager`](https://dumux.org/docs/doxygen/master/class_dumux_1_1_cake_grid_manager.html), which constructs a structured quadrilateral grid directly in polar coordinates: 80 radial cells per zone with mirrored grading toward both cylinder walls, and 320 uniform angular cells over the full $360°$ (`params.input`, used by the regression test). For the comparison, a second run with one global refinement of this grid (`-Grid.Refinement 1`, 320 × 640 cells) is performed.

**Result**

To run the benchmark and produce the plot and table below, execute:
```bash
python3 compile_run_plot.py
```
With `--update-readme`, the error table below and `images/analytical_comparison.png` are updated automatically.
The script builds the test, runs it on two grids and compares both with the analytical solution.
It expects PyVista and Matplotlib to be available for post-processing.
The coarse grid is the one used by the regression test (`params.input`), the finer grid is the same grid with one global refinement (`-Grid.Refinement 1`, every cell is split into four).

The script produces two files in the build directory of the test:
- `analytical_comparison.png`: radial velocity and pressure profiles of both numerical solutions and the analytical solution
- `l2_errors.md`: relative L2 errors of both runs

![Analytical comparison](images/analytical_comparison.png)

<!-- L2-ERRORS-START -->
| Grid (radial × angular cells) | Total cells | Rel. L2 error pressure | Rel. L2 error velocity |
|:--|--:|--:|--:|
| 160 × 320 | 51200 | 2.111e-03 | 9.079e-03 |
| 320 × 640 | 204800 | (copy from l2_errors.md) | (copy from l2_errors.md) |
<!-- L2-ERRORS-END -->
