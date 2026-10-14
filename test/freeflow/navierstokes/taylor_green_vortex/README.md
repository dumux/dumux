# Taylor-Green vortex {#benchmark-taylor-green-vortex}

This benchmark compares the spatial and temporal discretization schemes for the incompressible
Navier-Stokes equations available in DuMu<sup>x</sup> against an exact solution, the two-dimensional
Taylor-Green vortex. It is spatially periodic and decays exponentially in time.

## Problem description

### Governing equations

The incompressible Navier-Stokes equations read

$$
\begin{align*}
\rho \frac{\partial \mathbf{u}}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u})
- \nabla \cdot (\mu \nabla \mathbf{u}) + \nabla p &= \mathbf{f}, \\
\nabla \cdot \mathbf{u} &= 0,
\end{align*}
$$

| Symbol | Meaning |
|--------|---------|
| $\mathbf{u} = (u, v[, w])$ | velocity |
| $p$ | pressure |
| $\rho$ | density |
| $\mu = \rho \nu$, $\nu$ | dynamic and kinematic viscosity |
| $\mathbf{f}$ | source term (zero for the instationary Navier-Stokes problem) |

The viscous term is discretized in the form $\nabla \cdot (\mu \nabla \mathbf{u})$
(`FreeFlow.EnableUnsymmetrizedVelocityGradient = true`). For $\nabla \cdot \mathbf{u} = 0$ and constant $\mu$,
it equals the symmetrized form $\nabla \cdot (\mu (\nabla \mathbf{u} + \nabla \mathbf{u}^T))$, which is used with
`FreeFlow.EnableUnsymmetrizedVelocityGradient = false`.

### Construction of the solution

The solution belongs to the class of solutions described by Ethier & Steinman @cite EthierSteinman1994,
summarized in Antuono @cite Antuono2020 (Section 2):

* the velocity is divergence free,
* the time derivative balances the viscous term, $\rho \partial_t \mathbf{u} = \nabla \cdot (\mu \nabla \mathbf{u})$,
  so the solution decays exponentially in time,
* the convective term is a gradient (generalized Beltrami flow), $\nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p$,
  so it is balanced by the pressure alone.

Hence, the solution is exact for any amplitude and viscosity. It is periodic in every direction with the
wavelength $L$, i.e. with the wave number $k = 2\pi / L$.

### Two-dimensional case: Taylor-Green vortex

The bi-periodic solution of Taylor @cite Taylor1923Decay (in the form of Antuono @cite Antuono2020, Section 2.1,
scaled with the amplitude $U_0$) is an array of counter-rotating vortices:

$$
\begin{align*}
u(x, y, t) &= -U_0 \cos(k x) \sin(k y) F(t), \\
v(x, y, t) &= U_0 \sin(k x) \cos(k y) F(t), \\
p(x, y, t) &= -\frac{\rho U_0^2}{4} \left( \cos(2 k x) + \cos(2 k y) \right) F^2(t),
\end{align*}
$$

* decay factor: $F(t) = e^{-2 \nu k^2 t}$,
* vorticity (a scalar): $\omega = \partial_x v - \partial_y u = 2 k U_0 \cos(k x) \cos(k y) F(t)$.

This solution is due to Taylor (1923) alone. Taylor & Green @cite Taylor1937 studied the evolution of a
three-dimensional initial condition, for which no closed-form solution is known.

### Variants and source terms

The instationary Navier-Stokes problem uses **no source term**. The other variants omit a term of the
Navier-Stokes equations, which is therefore added as source term, evaluated with the analytical solution. With the
two identities above, all source terms are closed-form expressions in the analytical velocity and pressure:

| Variant | Navier-Stokes | Stokes (`Problem.EnableInertiaTerms = false`) |
|---------|---------------|-----------------------------------------------|
| instationary (`Problem.IsStationary = false`) | $\mathbf{f} = \mathbf{0}$ | $\mathbf{f} = -\nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = \nabla p$ |
| stationary (`Problem.IsStationary = true`) | $\mathbf{f} = -\rho \partial_t \mathbf{u} = \lambda \mu \mathbf{u}$ | $\mathbf{f} = \lambda \mu \mathbf{u} + \nabla p$ |

* $\lambda = 2 k^2$.
* The stationary variant is the solution at $t = 0$ ($F \equiv 1$): it adds the **storage term**.
* The Stokes variant adds the **convective term**.

### Initial and boundary conditions

Apart from the source terms, the analytical solution is only used for the initial and boundary conditions and the
errors:

* **Initial condition**: the $L^2$ projection of the analytical solution at $t = 0$ onto the discrete spaces. It is
  also the initial guess of the Newton solver for the stationary problem.
* **Velocity**: the analytical velocity is prescribed as Dirichlet condition on the entire boundary, evaluated at the
  current time (and at the stage times of the multi-stage methods).
* **Pressure**: since the velocity is prescribed on the entire boundary, the pressure is only determined up to a
  constant. It is therefore fixed to its analytical value at a single degree of freedom (the Box control volume of
  the vertex at the lower left corner).

Due to the periodicity of the solution, periodic boundary conditions on a domain of an integer number of
wavelengths would be the natural choice. They are not yet supported for the control-volume finite-element schemes
in the coupled free-flow setup.

## Discretization

### Momentum and mass balance

The mass balance is always discretized with the Box scheme (@ref BoxDiscretization). The momentum balance is
discretized with the following control-volume finite-element schemes:

| Scheme | Name | Grid | Remarks |
|--------|------|------|---------|
| PQ1Bubble (@ref PQ1BubbleDiscretization) | `pq1bubble` | simplices (`Dune::ALUGrid`) | one bubble per element, with the Box pressure only inf-sup stable on simplices |
| hybrid PQ1Bubble | `pq1bubblehybrid` | cubes (`Dune::YaspGrid`) | |
| hybrid PQ2 | `pq2hybrid` | cubes | |
| hybrid PQ2, higher-order quadrature | `pq2hybridhoquad` | cubes | see below |
| hybrid PQ3 | `pq3hybrid` | cubes | |

* In the hybrid schemes, only the vertex degrees of freedom have control volumes (the sub-control volumes of the
  Box scheme). The balance equations of the edge, face and element degrees of freedom are discretized with finite
  elements.
* The new CVFE problem interface is used (see also the tests `donea` and `sincos`).

### Quadrature rules

* The integrals are evaluated with the default quadrature rules of the schemes:
  * hybrid PQ2: midpoint rule on the sub-control volumes (only second-order consistent for a non-polynomial
    source), Gauss rule of order 2 on the sub-control-volume faces, order 4 for the finite-element terms;
  * hybrid PQ3: Gauss rules of order 3 for the sub-control volumes and faces, order 6 for the elements.
* `pq2hybridhoquad` checks whether the default rules of hybrid PQ2 are accurate enough. It uses a Gauss rule of
  order 3 on the sub-control volumes (`PQ2QuadratureTraits`) and of order 2 for the mass fluxes over the
  sub-control-volume faces of the Box scheme (`BoxQuadratureTraits`).
* The storage terms are lumped in all schemes (independent of the quadrature rules), which limits the spatial order
  of the instationary problem to 2.
* The errors and indicators are computed with a quadrature rule of order 8. For order 5, the quadrature error of
  the error norm of the higher-order schemes is of the same order as the error itself.

### Convective flux

* The convective flux over a sub-control-volume face transports the weighted two-point velocity
  $\zeta \mathbf{u}_\mathrm{up} + (1 - \zeta) \mathbf{u}_\mathrm{down}$ of the two adjacent sub-control volumes
  (`advectiveMomentumFlux` in `dumux/freeflow/navierstokes/momentum/cvfe/flux.hh`).
* The upwind weight is $\zeta$ = `Flux.UpwindWeight` $= 0.5$ (central, two-point average) in all studies and tests.
  Full upwinding ($\zeta = 1$) is only first-order accurate and adds numerical diffusion that visibly damps the
  vortex.
* The two-point average does not use the finite-element interpolant of the velocity at the face, which may limit
  the convergence order of the higher-order schemes. Possible remedies (not implemented yet):
  * evaluate the convected velocity with the finite-element interpolant at the integration points of the faces
    (as in the central scheme of the CVFEM in Nalu @cite Domino2018),
  * linear upwinding with consistently projected nodal gradients,
  * a deferred correction approach.

### Time integration

* The instationary problem is integrated with the multi-stage methods implicit Euler (default), Crank-Nicolson,
  DIRK2 and DIRK3 (`TimeLoop.Scheme`, see @ref benchmark-timestepping-methods) with a constant time step size.
* The multi-stage multi-domain assembler does not implement the assembly of the residual alone, so the Newton
  solver only uses the shift criterion (`Newton.EnableAbsoluteResidualCriterion = false`,
  `Newton.SatisfyResidualAndShiftCriterion = false`).
* The stationary problem is solved with the standard multi-domain assembler.

## Parameters

| Parameter | Symbol | Value | Unit | Input parameter |
|-----------|--------|-------|------|-----------------|
| Domain | $\Omega$ | $[0, 1]^2$ | m | `Grid.UpperRight` |
| Grid | | $32 \times 32$ | - | `Grid.Cells` |
| Wavelength | $L$ | $1$ | m | `Problem.WaveLength` (default: domain length in $x$ direction) |
| Wave number | $k = 2\pi/L$ | $2\pi$ | 1/m | |
| Amplitude | $U_0$ | $1$ | m/s | `Problem.ReferenceVelocity` |
| Root mean square velocity | $U = \left( \frac{1}{\|\Omega\|} \int_\Omega \|\mathbf{u}(\mathbf{x}, 0)\|^2 \, \mathrm{d}x \right)^{1/2}$ | $U_0/\sqrt{2} \approx 0.71$ | m/s | |
| Density | $\rho$ | $1$ | kg/m³ | `Component.LiquidDensity` |
| Kinematic viscosity | $\nu$ | $0.01$ | m²/s | `Component.LiquidKinematicViscosity` |
| Reynolds number | $Re = U L / \nu$ | $71$ | - | |
| End time | $T$ | $1$ | s | `TimeLoop.TEnd` |
| Time step size | $\Delta t$ | $0.05$ | s | `TimeLoop.DtInitial` |
| Time integration | | implicit Euler | - | `TimeLoop.Scheme` |
| Decay factor at end time | $F(T)$ | $0.45$ | - | |

* The parameters are set in `params_2d.input`.

## Benchmark indicators

All indicators are evaluated by quadrature (order 8) for the numerical solution $\mathbf{u}_h$, $p_h$ and the
analytical solution, and written for every time step to the error file (`Problem.ErrorFile`, default
`<Problem.Name>_errors.csv`).

### Vorticity, kinetic energy and enstrophy

| Indicator | Definition |
|-----------|------------|
| Vorticity | $\omega = \partial_x v - \partial_y u$ (a scalar) |
| Kinetic energy | $E(t) = \frac{1}{2} \int_\Omega \rho \|\mathbf{u}\|^2 \, \mathrm{d}x$ |
| Enstrophy | $\boldsymbol{\varepsilon}(t) = \frac{1}{2} \int_\Omega \rho \|\boldsymbol{\omega}\|^2 \, \mathrm{d}x$ |

* The vorticity of the numerical solution is evaluated with the gradient of the finite-element interpolant
  $\mathbf{u}_h$ (discontinuous across elements). It is also written to the ParaView output: at the element centers
  in the standard output and in the higher-order output.
* For an incompressible flow with periodic (or no-slip) boundaries, the kinetic energy is dissipated by viscosity at
  the rate $-\frac{\mathrm{d}E}{\mathrm{d}t} = 2 \nu \boldsymbol{\varepsilon}$.
* If $\Omega$ consists of an integer number of wavelengths in each direction (one wavelength by default), the
  analytical values are

  | | $E(t)$ | $\boldsymbol{\varepsilon}(t)$ |
  |-|--------|-------------------------------|
  | 2D | $\frac{1}{4} \rho U_0^2 \|\Omega\| F^2(t)$ | $2 k^2 E(t)$ |

  and the analytical solution satisfies the energy balance $-\frac{\mathrm{d}E}{\mathrm{d}t} = 2 \nu \boldsymbol{\varepsilon}$.
* For the numerical solution, the Dirichlet boundaries add boundary terms (the energy flux through $\partial \Omega$)
  and the time integration adds numerical dissipation.
* The solution, the identities behind the source terms and these formulas are verified symbolically by `verify.py`.

### Spatial convergence (stationary variant)

* Errors of velocity and pressure in the $L^2$ norm and the (full) $H^1$ norm,
  $\|\mathbf{u} - \mathbf{u}_h\|_{L^2(\Omega)}$, $\|\mathbf{u} - \mathbf{u}_h\|_{H^1(\Omega)}$,
  $\|p - p_h\|_{L^2(\Omega)}$, $\|p - p_h\|_{H^1(\Omega)}$, with
  $\|v\|_{H^1(\Omega)}^2 = \|v\|_{L^2(\Omega)}^2 + \|\nabla v\|_{L^2(\Omega)}^2$.
* The experimental orders of convergence (EOC) under uniform grid refinement are computed, as in
  @cite Schneider2024, with respect to the mesh sizes of velocity and pressure,
  $$
  h^v = \left( \frac{|\Omega|}{N_v} \right)^{1/d}, \qquad h^p = \left( \frac{|\Omega|}{N_p} \right)^{1/d}.
  $$
  $N_v$ and $N_p$ are the numbers of velocity and pressure degrees of freedom (nodes, not counting the velocity
  components; including bubble, edge and element degrees of freedom). This allows comparing schemes with different
  element types and numbers of degrees of freedom per element.

Expected orders:

| Momentum scheme | $\|\mathbf{u} - \mathbf{u}_h\|_{L^2(\Omega)}$ | $\|\mathbf{u} - \mathbf{u}_h\|_{H^1(\Omega)}$ | $\|p - p_h\|_{L^2(\Omega)}$ (Box) | $\|p - p_h\|_{H^1(\Omega)}$ (Box) |
|-----------------|---|---|---|---|
| PQ1Bubble | 2 | 1 | 2 | 1 |
| hybrid PQ1Bubble | 2 | 1 | 2 | 1 |
| hybrid PQ2 (default quadrature) | 2 | 2 | 2 | 1 |
| hybrid PQ2 (higher-order quadrature) | 2 | 2 | 2 | 1 |
| hybrid PQ3 | 3 | 2 | 2 | 1 |

* The Box pressure (P1/Q1) only converges with order 2, so the velocity of the hybrid schemes converges at most
  with order 2 in $H^1$ and 3 in $L^2$ (the schemes are not pressure-robust). Hybrid PQ3 is therefore not expected
  to reach the order 4 of a cubic interpolation (cf. the remark on P3-P1 in
  `test/freeflow/navierstokes/channel/3d_nonuniform`).
* The reduced velocity $L^2$ order 2 of hybrid PQ2 is a property of the hybrid P2 scheme itself (according to the
  developers of the hybrid schemes, due to cancellations in the integration), which hybrid PQ3 does not have.

### Temporal convergence (instationary variant)

* The errors at $t = T$ on a fine grid under refinement of the time step size.
* Expected orders: 1 (implicit Euler), 2 (Crank-Nicolson, DIRK2), 3 (DIRK3).

## How to reproduce results

### Executables

The executables are named `test_ff_navierstokes_taylorgreen_<dim>_<scheme>`:

* `<dim>`: `2d`,
* `<scheme>`: `pq1bubble`, `pq1bubblehybrid`, `pq2hybrid`, `pq2hybridhoquad` (hybrid PQ2 with higher-order quadrature
  rules) or `pq3hybrid`.

### Workflow

The benchmark consists of four scripts, which are executed in the build directory of this test
(`build-cmake/test/freeflow/navierstokes/taylor_green_vortex`):

```bash
cmake --build ../../../.. --target build_taylorgreen_tests
python3 verify.py                   # symbolic check of the analytical solution (requires sympy)
python3 benchmark.py --study all    # run the spatial, temporal and energy studies
python3 visualize.py                # plot the results (requires matplotlib, pyvista)
python3 benchmark.py --study animation --dim 2 && python3 animate.py  # animation (requires pyvista, ffmpeg)
```

| Script | Purpose | Output |
|--------|---------|--------|
| `verify.py` | symbolic check (`sympy`) of the solution, the source-term identities and the indicator formulas | console |
| `benchmark.py` | runs the studies | `runs/` |
| `visualize.py` | plots of all studies (matplotlib) | `results/taylorgreen_<study>_<dim>d.png` |
| `animate.py` | animation of the fields of `benchmark.py --study animation`, sampled from the higher-order output with [pyvista](https://pyvista.org), plotted with matplotlib and written with ffmpeg | `results/taylorgreen_animation_2d_<scheme>.mp4` |

### Output of `benchmark.py`

| Path | Content |
|------|---------|
| `runs/<name>/<name>_errors.csv` | errors and indicators of every time step of a run |
| `runs/<name>/vtk/<name>.pvd` | standard VTK output (Box pressure, element velocities and vorticity, analytical solution) |
| `runs/<name>/vtk/<name>_ho.pvd` | higher-order output: velocity, pressure and vorticity on Lagrange cells of the polynomial degree of the velocity (requires GridFormat; increase the "Nonlinear Subdivision Level" in ParaView to render the curved fields) |
| `runs/taylorgreen_<study>_<dim>d.md` | convergence tables |
| `runs/taylorgreen_<study>_<dim>d.json` | manifest of the runs of a study for `visualize.py` and `animate.py` |

By default, ParaView output is written for every time step of the instationary runs (temporal and energy study,
animation).

### Options of `benchmark.py`

| Option | Meaning | Default |
|--------|---------|---------|
| `--study spatial temporal energy animation all` | studies to run (`all`: all except the animation) | `all` |
| `--schemes ...` | momentum schemes | spatial and energy: all; temporal: hybrid PQ2; animation: hybrid PQ3 |
| `--time-scheme ImplicitEuler\|CrankNicolson\|DIRK2\|DIRK3` | time integration of the energy study and the animation | `ImplicitEuler` |
| `--wavelength L` | wavelength of the solution | from the input file |
| `--upwind-weight W` | upwind weight of the convective flux (`Flux.UpwindWeight`), e.g. `1` for full upwinding | 0.5 |
| `--cells N1 N2 ...` | cells per direction of the grids of the spatial study | 16, 32, 64, 128 (PQ1Bubble: 32, 64, 128) |
| `--levels N` | number of time step refinements of the temporal study | 4 on $64^2$ cells |
| `--vtk` | ParaView output for all runs | instationary runs only |
| `--reuse` | reuse existing error files instead of rerunning the simulations | |
| `--test` | short study checking the convergence rates (used by CTest) | |

### Options of `visualize.py` and `animate.py`

| Option | Meaning |
|--------|---------|
| `--study spatial temporal energy` | plots to create (`visualize.py`) |
| `--output-dir` | directory of the figures and the video (default: `results`) |
| `--resolution`, `--arrows`, `--fps` | sampling points and arrows per direction, frames per second (`animate.py`) |

## Related setups and literature

### Related setups

* The stationary three-dimensional Stokes test of Schneider & Koch @cite Schneider2024 uses the Taylor-Green initial
  condition $\mathbf{u} = (-2\cos\sin\sin, \sin\cos\sin, \sin\sin\cos)$ together with a manufactured pressure and
  source. The convergence plots of this benchmark follow their definition of the mesh sizes $h^v$ and $h^p$.
* The [OpenFOAM verification and validation case](https://gitlab.com/openfoam/community/verification-and-validation/-/tree/main/cases/Taylor-Green-vortex)
  simulates the two-dimensional vortex on a periodic domain $[0, 2\pi]^2$ with $k = 1$, $\nu = 1$, $\Delta t = 0.05$
  until $t = 0.3$ and compares velocity and pressure profiles. An equivalent setup (up to the boundary conditions and
  the sign of the velocity) is obtained with
  `-Grid.UpperRight "6.283185307179586 6.283185307179586" -Component.LiquidKinematicViscosity 1 -TimeLoop.TEnd 0.3`.

### Schemes used by other groups

Combinations of schemes used for Taylor-Green(-type) verification:

| Reference | Velocity / pressure | Convective term | Time integration | Setup | Reported orders |
|-----------|---------------------|-----------------|------------------|-------|-----------------|
| @cite Kim1985 | staggered (MAC) finite differences | central, explicit (Adams-Bashforth 2) | fractional step, Crank-Nicolson for viscous term | 2D, instationary | 2 (space and time) |
| @cite Mortensen2015 (Oasis, FEniCS) | Taylor-Hood P2/P1, equal order P1/P1 | Adams-Bashforth extrapolated convecting velocity | incremental pressure correction (IPCS), Crank-Nicolson | 2D periodic, $\nu = 0.01$ | P2/P1: velocity 4 (superconvergence on uniform mesh, 3 otherwise), pressure 2; P1/P1: 2/2; time 2 |
| @cite Fehn2017, @cite Fehn2018 | discontinuous Galerkin $(k, k-1)$ | local Lax-Friedrichs flux | BDF1/2 splitting and coupled schemes | 2D vortex (Stokes and Navier-Stokes); 3D $Re = 1600$ (energy and dissipation only) | velocity $k+1$, pressure $k$ |
| @cite Domino2018 (Nalu) | vertex-centered CVFEM P1, P2 (equal order, pressure stabilization) | central with finite-element interpolation at faces | backward Euler, BDF2 | 2D convecting Taylor vortex | P1: 2, P2: 3 (space); 1, 2 (time) |
| Schneider & Koch @cite Schneider2024 | CVFEM P1-bubble/P1 (MINI-type, also hybrid) | none (Stokes) | stationary | 3D, Stokes | velocity $L^2$ 2 (2.2-2.5), $H^1$ 1, pressure $L^2$ 1.5 (Delaunay grids) |
| @cite Schneider2023 | staggered finite volumes (coupled with Box) | weighted upwinding $\zeta$ | implicit Euler | 2D manufactured, coupled free flow/porous medium | about 2 (velocity and pressure) |

### Pressure robustness

* The convective term of the solution is a gradient balanced by the pressure.
* For pairs that are not pressure-robust (e.g. MINI, Taylor-Hood, and also the schemes with Box pressure used here),
  the velocity error then contains a contribution proportional to the pressure approximation error divided by $\nu$
  @cite John2017, @cite Gauger2019, which grows with the Reynolds number. The convergence studies therefore use moderate
  Reynolds numbers.
* Alternative formulations of the convective term with improved conservation properties (e.g. EMAC) are discussed
  in @cite Charnyi2017.
