# Decaying Vortices (FVCA8 Benchmark) {#benchmark-angeli}

## Unsteady Navier-Stokes flow with an analytical solution

**Problem Description**

This benchmark is the unsteady two-dimensional test of the FVCA8 benchmark for the Stokes and Navier-Stokes
equations @cite Boyer2016, for which Angeli et al. @cite Angeli2017 report results. It solves the unsteady,
incompressible Navier-Stokes equations without source terms,

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = -\nabla p + \nu \nabla^2 \mathbf{u},
\qquad \nabla \cdot \mathbf{u} = 0,
\qquad \int_\Omega p \, \mathrm{d}\mathbf{x} = 0,
$$

for the velocity $\mathbf{u} = (u_x, u_y)$ and the pressure $p$ on the unit square, with the kinematic
viscosity $\nu$ and the density $\rho = 1$ kg/m³. The solution is a pattern of vortices whose shape does not
change while it decays exponentially in time,

$$
u_x = -2\pi\, e^{-5\nu\pi^2 t} \cos(\pi x) \sin(2\pi y), \qquad
u_y = \pi\, e^{-5\nu\pi^2 t} \sin(\pi x) \cos(2\pi y),
$$

$$
p = -\frac{\pi^2}{4}\, e^{-10\nu\pi^2 t} \left( 4\cos(2\pi x) + \cos(4\pi y) \right).
$$

The advective term does not vanish, so the benchmark also tests its discretization. Since no source term is
needed, the equations are solved exactly as they are used in applications.

**Parameters**

| Parameter           | Symbol   | Value  | Unit  |
|---------------------|----------|--------|-------|
| Domain              | $\Omega$ | $(0, 1)^2$ | m |
| Fluid density       | $\rho$   | $1$    | kg/m³ |
| Kinematic viscosity | $\nu$    | $0.1$, $0.01$ | m²/s |
| Final time          | $T = 1/(10\nu)$ | $1$, $10$ | s |

The largest initial velocity $2\pi$ m/s corresponds to Reynolds numbers of about 63 and 628, based on the side
length of the domain. At the final time, the velocity has decayed to 0.7 % and the pressure to 0.005 % of its
initial value.

**Setup**

The velocity is set to the analytical solution on the entire boundary, and the initial velocity is the analytical
solution, as prescribed by the benchmark. The benchmark leaves the choice of the scheme to the participants:

- **Space:** The momentum balance is discretized with the face-centered staggered scheme
  (@ref FaceCenteredStaggeredDiscretization) and the mass balance with the cell-centered two-point flux
  approximation (@ref CCTpfaDiscretization), on the uniform Cartesian grids `mesh_cart_1` to `mesh_cart_5` of the
  benchmark. With $n$ cells per direction, there are $2n(n + 1)$ velocity and $n^2$ pressure unknowns. The finer
  grids `mesh_cart_6` and `mesh_cart_7` are not used because of their runtime.
- **Advection:** The advective momentum flux uses central weighting (`Flux.UpwindWeight = 0.5`), or full upwinding
  (`Flux.UpwindWeight = 1`) for comparison.
- **Time:** The implicit Euler method is used with the fixed time step size $\Delta t = 0.2/(\nu n^2)$. It decreases
  with the square of the cell size, like the expected spatial error.
- **Nonlinear solver:** The Newton method stops when the maximum relative shift is below $10^{-5}$.
  The tables list the mean number of iterations per time step.
- **Pressure:** The pressure is fixed to the analytical solution in one cell and made mean-free before the errors
  are computed.
- **Velocity data:** The initial and boundary velocities are averages of the analytical solution over the faces.

These settings are passed on the command line. The regression tests do not use them and rely on the settings in
`params.input`, with the default full upwinding, a short simulation time and the pressure fixed in all cells that
touch the boundary.

**Result**

The figure and the tables are produced by running
```bash
python3 plot_results.py ./test_ff_navierstokes_angeli
```
in the build directory `test/freeflow/navierstokes/angeli` (requires NumPy, Matplotlib and PyVista). The
script takes about 40 minutes on eight cores.

![Velocity and pressure on a vertical line at the final time for both weightings (top), and the errors for both viscosities under refinement (bottom)](angeli_comparison.png){html: width=90%}

The tables list the relative space-time $L^2$ errors and the orders defined by the benchmark,

$$
\mathrm{erru} = \left[ \frac{\int_0^T \int_\Omega |\mathbf{u} - \mathbf{u}_\mathrm{ex}|^2}{\int_0^T \int_\Omega |\mathbf{u}_\mathrm{ex}|^2} \right]^{1/2},
\qquad
\mathrm{ordu}_i = -2\, \frac{\ln \mathrm{erru}_i - \ln \mathrm{erru}_{i-1}}{\ln \mathrm{nuu}_i - \ln \mathrm{nuu}_{i-1}},
$$

and analogously for the pressure with the number of pressure unknowns and for the velocity gradient (`errgu`). Here,
$\mathrm{nuu}$ is the number of velocity unknowns. The integrals are evaluated at the degrees of freedom and the time
steps. The velocity gradient is evaluated with the difference quotients of the viscous flux, compared with the
analytical gradient at the midpoint of the two positions. The divergence error is the absolute norm
$\mathrm{errdivu} = \left[ \int_0^T \int_\Omega |\nabla \cdot \mathbf{u}|^2 \right]^{1/2}$ of the discrete divergence in
the cells. It stays at round-off level (at most about $3 \cdot 10^{-13}$) because the scheme satisfies the discrete mass
balance, so its order is not meaningful and is omitted.

For $\nu = 0.1$ m²/s and central weighting:

| Mesh | Cells | errgu | ordgu | erru | ordu | errp | ordp | errdivu | Newton iterations |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 4 × 4 | 1.721e-01 | - | 1.254e-01 | - | 4.160e-01 | - | 4.2e-16 | 3.1 |
| 2 | 8 × 8 | 6.851e-02 | 1.44 | 3.792e-02 | 1.87 | 1.080e-01 | 1.95 | 1.1e-15 | 3.0 |
| 3 | 16 × 16 | 1.807e-02 | 2.01 | 9.245e-03 | 2.12 | 2.771e-02 | 1.96 | 1.5e-15 | 2.8 |
| 4 | 32 × 32 | 4.538e-03 | 2.04 | 2.297e-03 | 2.05 | 6.991e-03 | 1.99 | 2.9e-15 | 2.5 |
| 5 | 64 × 64 | 1.136e-03 | 2.02 | 5.736e-04 | 2.02 | 1.752e-03 | 2.00 | 5.4e-15 | 2.3 |

For $\nu = 0.1$ m²/s and full upwinding:

| Mesh | Cells | errgu | ordgu | erru | ordu | errp | ordp | errdivu | Newton iterations |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 4 × 4 | 1.722e-01 | - | 1.103e-01 | - | 3.685e-01 | - | 2.6e-15 | 3.5 |
| 2 | 8 × 8 | 1.364e-01 | 0.36 | 7.519e-02 | 0.60 | 1.660e-01 | 1.15 | 2.8e-15 | 3.1 |
| 3 | 16 × 16 | 9.639e-02 | 0.52 | 3.962e-02 | 0.96 | 1.064e-01 | 0.64 | 1.5e-15 | 2.8 |
| 4 | 32 × 32 | 4.916e-02 | 0.99 | 1.859e-02 | 1.12 | 6.685e-02 | 0.67 | 3.0e-15 | 2.6 |
| 5 | 64 × 64 | 2.262e-02 | 1.13 | 9.258e-03 | 1.02 | 3.745e-02 | 0.84 | 5.3e-15 | 2.3 |

For $\nu = 0.01$ m²/s and central weighting:

| Mesh | Cells | errgu | ordgu | erru | ordu | errp | ordp | errdivu | Newton iterations |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 4 × 4 | 1.304e+00 | - | 1.036e+00 | - | 7.599e-01 | - | 1.1e-14 | 3.8 |
| 2 | 8 × 8 | 3.699e+01 | -5.22 | 1.539e+01 | -4.21 | 3.137e+02 | -8.69 | 3.1e-13 | 5.6 |
| 3 | 16 × 16 | 4.636e-01 | 6.59 | 1.161e-01 | 7.35 | 7.147e-02 | 12.10 | 7.8e-14 | 2.8 |
| 4 | 32 × 32 | 4.543e-02 | 3.43 | 8.243e-03 | 3.90 | 9.998e-03 | 2.84 | 9.2e-15 | 2.6 |
| 5 | 64 × 64 | 7.919e-03 | 2.55 | 1.647e-03 | 2.35 | 2.477e-03 | 2.01 | 1.7e-14 | 2.3 |

For $\nu = 0.01$ m²/s and full upwinding:

| Mesh | Cells | errgu | ordgu | erru | ordu | errp | ordp | errdivu | Newton iterations |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 4 × 4 | 3.293e-01 | - | 2.484e-01 | - | 4.245e-01 | - | 2.2e-14 | 3.5 |
| 2 | 8 × 8 | 3.899e-01 | -0.26 | 2.118e-01 | 0.25 | 2.566e-01 | 0.73 | 2.3e-14 | 3.1 |
| 3 | 16 × 16 | 4.282e-01 | -0.14 | 1.348e-01 | 0.68 | 1.526e-01 | 0.75 | 7.4e-15 | 2.8 |
| 4 | 32 × 32 | 3.376e-01 | 0.35 | 6.520e-02 | 1.07 | 8.580e-02 | 0.83 | 2.0e-14 | 2.6 |
| 5 | 64 × 64 | 1.865e-01 | 0.87 | 2.791e-02 | 1.24 | 4.815e-02 | 0.83 | 1.7e-14 | 2.3 |

The numbers of unknowns and of entries in the blocks of the Jacobian (the velocity-velocity block `nnzu`, the
pressure-pressure block `nnzp` and the velocity-pressure block `nnzup`) are the same for all cases. These are the
allocated entries of the sparsity pattern, some of which can be numerically zero:

| Mesh | Cells | nuu | npu | nnzu | nnzp | nnzup |
|--:|:--|--:|--:|--:|--:|--:|
| 1 | 4 × 4 | 40 | 16 | 292 | 64 | 64 |
| 2 | 8 × 8 | 144 | 64 | 1164 | 288 | 256 |
| 3 | 16 × 16 | 544 | 256 | 4636 | 1216 | 1024 |
| 4 | 32 × 32 | 2112 | 1024 | 18492 | 4992 | 4096 |
| 5 | 64 × 64 | 8320 | 4096 | 73852 | 20224 | 16384 |

Central weighting converges with second order, also in the velocity gradient. For $\nu = 0.01$ m²/s, it does so only on the grids with 32 × 32
cells and more: on the coarsest grids, the cell Péclet number is large and the solution is not resolved, with
errors above one and, on the grid with 8 × 8 cells, Newton failures that reduce the time step size. The orders on
the following grids are correspondingly too large. Full upwinding is stable on all grids, but converges with about
first order in velocity and pressure. On the finest grid, its errors in velocity and pressure are 16 to 21 times
larger than with central weighting, and its gradient errors are 20 to 24 times larger.

The pressure in the figure deviates visibly from the analytical solution at the final time, for both weightings.
The analytical pressure decays with the square of the velocity, whereas the errors of the linear terms in the
momentum balance decay only with the velocity, so the relative pressure error grows with time. For
$\nu = 0.1$ m²/s and central weighting, the relative pressure error at the final time is 1.95, 0.49 and 0.12 on the
grids with 16, 32 and 64 cells per direction: large, but of second order. The space-time errors in the tables weight
the early times, when the solution is largest.

The averaged velocity data is important for the pressure. Pointwise values of the analytical solution are not
discretely divergence-free, which causes a spurious pressure in the first time step. For $\nu = 0.1$ m²/s and central
weighting, this reduces `ordp` to 1.80, 1.52 and 1.25 on the meshes 2 to 4, while `ordu` remains about two.

**Comparison with Angeli et al.**

Angeli et al. @cite Angeli2017 report results for the same grids with the code TrioCFD, using a projection method,
an explicit third-order Runge-Kutta scheme and the QUICK scheme for advection on Cartesian grids. On the mesh 5:

| $\nu$ (m²/s) | Scheme | errgu | ordgu | erru | ordu | errp | ordp |
|--:|:--|--:|--:|--:|--:|--:|--:|
| 0.1 | Angeli et al. | 7.747e-03 | 1.29 | 2.148e-03 | 1.54 | 6.300e-03 | 1.52 |
| 0.1 | DuMux, central | 1.136e-03 | 2.02 | 5.736e-04 | 2.02 | 1.752e-03 | 2.00 |
| 0.1 | DuMux, upwind | 2.262e-02 | 1.13 | 9.258e-03 | 1.02 | 3.745e-02 | 0.84 |
| 0.01 | Angeli et al. | 1.362e-02 | 0.99 | 5.668e-03 | 1.56 | 7.711e-03 | 0.83 |
| 0.01 | DuMux, central | 7.919e-03 | 2.55 | 1.647e-03 | 2.35 | 2.477e-03 | 2.01 |
| 0.01 | DuMux, upwind | 1.865e-01 | 0.87 | 2.791e-02 | 1.24 | 4.815e-02 | 0.83 |

The errors of Angeli et al. are printed without the time integral of the benchmark definition, and their time step
size and discrete gradient are not specified, so the comparison is only indicative. Their numbers of velocity unknowns
are twice the ones given here, with the same ratios between the grids. Their `nnzp` equals the values above, for the
matrix of a pressure Poisson equation with the same five-point stencil.

**Tests**

The test `test_ff_navierstokes_angeli_convergence` runs the benchmark setup with central weighting and
$\nu = 0.1$ m²/s on the meshes 1 to 3, and checks that the mean orders of pressure and velocity are at least 1.8 and
that `errdivu` is below $10^{-10}$. The tests `test_ff_navierstokes_angeli` and `test_ff_navierstokes_angeli_averaged` compare the solution with reference
solutions.

**Variants**

With `-Problem.UseNeumann true`, the momentum flux of the analytical solution is prescribed on the right and top
boundaries, and the pressure is not fixed.
