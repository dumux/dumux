# Taylor-Couette Flow {#benchmark-taylor-couette}

## Steady viscous flow between two rotating cylinders

**Problem Description**

An incompressible fluid with density $\rho$ and dynamic viscosity $\mu$ fills the gap between two long,
concentric cylinders. The inner cylinder of radius $R_1$ rotates with the angular velocity $\Omega_1$,
the outer cylinder of radius $R_2$ with $\Omega_2$. If the flow does not vary along the axis, the steady
Navier-Stokes equations

$$
\rho (\mathbf{u} \cdot \nabla) \mathbf{u} = -\nabla p + \mu \nabla^2 \mathbf{u}, \qquad \nabla \cdot \mathbf{u} = 0
$$

for the velocity $\mathbf{u}$ and the pressure $p$ have the solution $u_r = 0$ and

$$
u_{\theta}(r) = A r + \frac{B}{r}, \qquad
A = \frac{\Omega_2 R_2^2 - \Omega_1 R_1^2}{R_2^2 - R_1^2},
\qquad
B = \frac{(\Omega_1 - \Omega_2) R_1^2 R_2^2}{R_2^2 - R_1^2},
$$

where $r$ is the distance from the axis and $u_r$, $u_\theta$ are the radial and tangential velocity.
The pressure follows from the radial momentum balance $\partial p/\partial r = \rho\, u_\theta^2/r$,

$$
p(r) = \rho \left( \frac{A^2 r^2}{2} + 2AB \ln r - \frac{B^2}{2r^2} \right) + C.
$$

The constant $C$ is arbitrary because the velocity is prescribed on the entire boundary.
Taylor @cite taylor1923 showed that this flow becomes unstable to axisymmetric vortices (Taylor vortices)
above a critical Reynolds number $Re = \Omega_1 R_1 (R_2 - R_1) \rho/\mu$. For $R_1/R_2 = 0.5$ and
$\Omega_2 = 0$, linear stability analysis against axisymmetric perturbations gives $Re_c = 68.2$.
The test solves the two-dimensional problem in the plane perpendicular to the axis. This model has no
axial direction and cannot represent Taylor vortices, so at $Re = 100$ it computes the circular Couette
flow above.

**Parameters**

| Parameter                       | Symbol     | Value  | Unit  |
|---------------------------------|------------|--------|-------|
| Inner cylinder radius           | $R_1$      | $1$    | m     |
| Outer cylinder radius           | $R_2$      | $2$    | m     |
| Inner cylinder angular velocity | $\Omega_1$ | $100$  | rad/s |
| Outer cylinder angular velocity | $\Omega_2$ | $0$    | rad/s |
| Fluid density                   | $\rho$     | $1$    | kg/m³ |
| Fluid dynamic viscosity         | $\mu$      | $1$    | Pa s  |

With these values, $Re = 100$.

**Setup**

The momentum balance is discretized with the face-centered diamond scheme (@ref DiamondDiscretization),
the mass balance with the cell-centered two-point flux approximation (@ref CCTpfaDiscretization).
The advective momentum flux uses central weighting (`Flux.UpwindWeight = 0.5`).
The annulus is meshed with the `CakeGridManager` as a structured grid of quadrilaterals with
16 cells over the gap and 128 cells in angular direction (`params.input`). Global refinement
(`-Grid.Refinement n`) halves the cell size in both directions. The new vertices lie on the edges of the
coarse grid, so the cylinders remain polygons with 128 edges.

On both cylinders, the velocity is set to the analytical solution, evaluated on the boundary of the grid.
Since the analytical solution also solves the problem in the polygonal domain covered by the grid,
the comparison measures the discretization error and not the approximation of the cylinders by polygons.
The pressure is fixed in one cell. After the solve, it is shifted such that its volume
average equals that of the analytical solution, so that the error in the fixed cell does not offset the
entire pressure field.

**Result**

The figure and the table are produced by running
```bash
python3 plot_results.py ./test_ff_navierstokes_taylorcouette
```
in the build directory `test/freeflow/navierstokes/taylorcouette` (requires NumPy, Matplotlib and PyVista).
The script runs the grid of `params.input` and three global refinements of it.

![Tangential velocity and pressure over the radius (top) and their errors (bottom) on all four grids](taylorcouette_comparison.png){html: width=90%}

The figure shows the values in the cells of one angular column. With each refinement, the error of both
velocity and pressure decreases by a factor of about four at every radius. The pressure error changes sign
near $r = 1.4$ m: after matching the volume-averaged pressure, its volume average is zero.

The table lists the discrete $L^2$ errors $\left( \sum_i |V_i|\, |e_i|^2 / \sum_i |V_i| \right)^{1/2}$,
where $e_i$ is the error at the degree of freedom $i$ with control volume $V_i$: the cell centers for the
pressure and the face centers for the velocity. For comparison, the analytical pressure difference across
the gap is $p(R_2) - p(R_1) = 2172$ Pa and the velocity of the inner cylinder is $\Omega_1 R_1 = 100$ m/s.

| Grid (radial × angular cells) | $L^2$ error pressure (Pa) | rate | $L^2$ error velocity (m/s) | rate |
|:--|--:|--:|--:|--:|
| 16 × 128 | 1.519e+02 | - | 6.727e+00 | - |
| 32 × 256 | 3.368e+01 | 2.17 | 1.534e+00 | 2.13 |
| 64 × 512 | 8.143e+00 | 2.05 | 3.736e-01 | 2.04 |
| 128 × 1024 | 2.018e+00 | 2.01 | 9.277e-02 | 2.01 |

Both errors converge with second order. The test `test_ff_navierstokes_taylorcouette` compares the solution
on the coarsest grid with a reference solution, and `test_ff_navierstokes_taylorcouette_convergence` checks that
the mean convergence rate over the first three grids is at least 1.9.
