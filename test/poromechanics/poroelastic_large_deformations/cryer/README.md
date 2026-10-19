# Cryer's Sphere Consolidation {#benchmark-cryer}

A fluid-saturated porous sphere is loaded on its drained surface. The pore pressure at the centre
first rises above the load and then decays (Mandel–Cryer effect). Cryer @cite Cryer1963 solved the
problem in linear poroelasticity. The benchmark solves it with linear poroelasticity
(@ref PoroElastic coupled to @ref OnePModel, `test/poromechanics/poroelastic/cryer`) and with the
large-deformation poroelastic model (@ref PoroElasticLargeDef,
`test/poromechanics/poroelastic_large_deformations/cryer`). With a small load, both are compared with
Cryer's solution. With larger loads, the large-deformation model shows the effect of finite deformation
and of a permeability that depends on the deformation.

## Problem

The sphere has the radius $R_0 = 1$ m; one octant is modelled.
- Sphere surface: the pressure load $p_0$, the traction $-p_0 \mathbf{N}$ with the radial unit vector
  $\mathbf{N}$ of the reference configuration, and the drained condition $p_f = 0$.
- Symmetry planes $x = 0$, $y = 0$, $z = 0$: zero normal displacement and no flow.

![Section z = 0 of the simulated octant: load and drainage on the sphere surface, rollers on the symmetry planes x = 0 and y = 0](cryer_setup.svg){html: width=50%}

In linear poroelasticity, the total stress is $\boldsymbol{\sigma} = \lambda \operatorname{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2G \boldsymbol{\varepsilon} - \alpha_B p_f \mathbf{I}$
with $\lambda = K - \tfrac{2}{3} G$, and the porosity $\phi = \phi_0 + \alpha_B \nabla \cdot \mathbf{u} + S_p p_f$ gives the
storage $\alpha_B \partial_t \nabla \cdot \mathbf{u} + S_p \partial_t p_f$.

In the large-deformation model, the skeleton is a Neo-Hookean material with volumetric–isochoric
split, isochoric stress $\mathbf{P}_\text{iso} = G J^{-2/3}(\mathbf{F} - \tfrac{1}{3}\operatorname{tr}(\mathbf{C})\, \mathbf{F}^{-\top})$ and
solid pressure $p_s = \frac{K}{2}(J - J^{-1})$. The model adds $(p_s - \alpha_B p_f) J \mathbf{F}^{-\top}$ to the
stress and has the storage $\alpha_B \partial_t J + S_p \partial_t p_f$; see @ref PoroElasticLargeDef for the
equations. For small deformations, both models have the same linearization.

| Parameter | Symbol | Value |
| :--- | :--- | ---: |
| shear modulus | $G$ | 1.5 MPa |
| drained bulk modulus | $K$ | 1.0 MPa |
| Poisson's ratio of the linearized skeleton | $\nu$ | 0 |
| Biot coefficient | $\alpha_B$ | 1 |
| storage coefficient | $S_p$ | 0 |
| initial porosity | $\phi_0$ | 0.5 |
| initial permeability | $\kappa_0$ | $10^{-12}$ m² |
| fluid viscosity | $\mu_f$ | $10^{-3}$ Pa s |
| surface load | $p_0$ | 100 Pa, 0.25 MPa or 0.5 MPa ($p_0/K = 10^{-4}$, 0.25, 0.5) |

The permeability of the large-deformation model is either constant, $\kappa = \kappa_0$, or follows the
Kozeny–Carman law $\kappa(J) = \frac{d_p^2}{180} \frac{\phi^3}{(1 - \phi)^2}$ with the porosity
$\phi(J) = 1 - (1 - \phi_0)/J$ and the particle diameter $d_p$ chosen such that $\kappa(1) = \kappa_0$
(`SpatialParams.PermeabilityModel`). In linear poroelasticity, the permeability is constant, and
$p_f/p_0$ does not depend on the load.

## Cryer's solution

In linear poroelasticity, the pore pressure is

$$\frac{p(r,t)}{p_0} = 2mB \sum_{i=1}^{\infty}
   \frac{\sin z_i - z_i \operatorname{sinc}(z_i r/R_0)}{m z_i \cos z_i + (2m-1)\sin z_i}
   \exp\left(-\frac{z_i^2 c_v t}{\beta R_0^2}\right), \qquad \operatorname{sinc}(x) = \frac{\sin x}{x},$$

where the $z_i$ are the positive roots of $(1 - m z^2)\sin z - z\cos z = 0$ and

$$m = \frac{K + \tfrac{4}{3}G}{4G}\left(1 + \frac{K S_p}{\alpha_B^2}\right), \qquad
c_v = \frac{\kappa}{\mu_f}\left(K + \tfrac{4}{3}G\right), \qquad
\beta = \alpha_B^2 + \left(K + \tfrac{4}{3}G\right) S_p, \qquad
B = \frac{\alpha_B}{\alpha_B^2 + K S_p}.$$

$B$ is Skempton's coefficient, the uniform pore pressure $B p_0$ right after loading. With the
parameters above, $m = 0.5$, $c_v = 3 \times 10^{-3}$ m²/s, $\beta = 1$, $B = 1$, and the centre
pressure has its maximum 1.575 $p_0$ at $t c_v/R_0^2 = 0.0625$. `test/poromechanics/cryer_analytical.hh` implements the series.

## Discretization

Linear poroelasticity: the pressure is cell-centred with the multi-point flux approximation
(@ref CCMpfaDiscretization), since the two-point flux approximation is not consistent on the
tetrahedral mesh; the displacement is linear with the box scheme (@ref BoxDiscretization). The two
subdomains are coupled with `PoroMechanicsCouplingManager` (@ref MultiDomain). Permeability and fluid
density do not depend on the deformation, so the coupling manager of the test restricts the dependence
of the flow residual on the displacement to the vertices of the element (`test/poromechanics/poroelastic/cryer/couplingmanager.hh`).
The centre pressure is the value of the cell at the centre, compared with Cryer's solution at the
centre of that cell.

Large deformation: the displacement is quadratic with the hybrid control-volume / finite-element
scheme (@ref PQ2Discretization), solid and fluid pressure are linear with the box scheme
(@ref BoxDiscretization), coupled as three subdomains (@ref MultiDomain).

The tests use the mesh `test/poromechanics/sphere_eighth_coarse.msh` (223 vertices, 670 tetrahedra). The
benchmark runs use the mesh `sphere_eighth.msh` that gmsh generates from `test/poromechanics/sphere_eighth.geo`
with the maximum element size
0.095 m; with gmsh 4.13.1 it has 862 vertices and 3297 tetrahedra. Time integration uses the multi-stage
time stepper with the third-order diagonally implicit Runge–Kutta method DIRK3 (`TimeLoop.Scheme`,
implicit Euler is also available). The time step starts at 1 s and grows by the factor 1.3 per step up
to 20 s, which reaches $R_0^2/c_v = 333.33$ s in 25 steps.

## Results

Centre pressure $p/p_0$ and volume ratio $J$ at the centre, mesh `sphere_eighth.msh` (gmsh 4.13.1), DIRK3,
25 steps:

| model | permeability | $p_0/K$ | maximum of $p/p_0$ (at $t c_v/R_0^2$) | $J$ at $t c_v/R_0^2 = 1$ | $p/p_0$ at $t c_v/R_0^2 = 1$ |
| :--- | :--- | ---: | ---: | ---: | ---: |
| linear | constant | any | 1.565 (0.072) | 0.9999 | 0.031 |
| large deformation | constant | $10^{-4}$ | 1.583 (0.072) | 0.9999 | 0.031 |
| large deformation | constant | 0.25 | 1.710 (0.072) | 0.743 | 0.028 |
| large deformation | constant | 0.5 | 1.865 (0.072) | 0.473 | 0.017 |
| large deformation | Kozeny–Carman | $10^{-4}$ | 1.583 (0.072) | 0.9999 | 0.031 |
| large deformation | Kozeny–Carman | 0.25 | 1.587 (0.072) | 0.802 | 0.505 |
| large deformation | Kozeny–Carman | 0.5 | 1.565 (0.096) | 0.738 | 1.058 |

With the small load, the centre pressure deviates most from Cryer's solution at the first time step,
$t c_v/R_0^2 = 0.003$: by $-0.026\,p_0$ in linear poroelasticity and by $0.018\,p_0$ in the large-deformation
model; after $t c_v/R_0^2 = 0.1$, the deviations are at most $0.0025\,p_0$ and $0.0075\,p_0$. At
$t c_v/R_0^2 = 0.072$, Cryer's solution is 1.568 $p_0$ at the centre and 1.562 $p_0$ at the centre of the
centre cell of linear poroelasticity. With the larger loads, the
overshoot grows for a constant permeability. With the Kozeny–Carman law, the compaction reduces the
permeability, for $J = 0.738$ to 15 % of $\kappa_0$, and the centre has not drained at
$t c_v/R_0^2 = 1$.

![Centre pressure (top) and volume ratio at the centre (bottom) for the three loads of the large-deformation model, constant (left) and Kozeny–Carman permeability (right); dash-dotted: linear poroelasticity, dashed: Cryer's solution](cryer_overview.png)

## Tests

- `test_cryer_poroelastic_small_def`: linear poroelasticity, coarse mesh, implicit Euler, five steps of
  1 s. The test fails if the centre pressure deviates from Cryer's solution by more than 0.07 $p_0$
  (`Problem.CenterPressureTolerance`); the deviation is 0.0602 $p_0$.
- `test_cryer_poroelastic_small_def_compressible`: as above with $\alpha_B = 0.8$ and
  $S_p = 5 \times 10^{-7}$ 1/Pa ($B = 0.70$), to 180 s with growing time steps; the deviation is
  0.0397 $p_0$.
- `test_cryer_poroelastic_large_def`: large-deformation model with the small load, otherwise as
  `test_cryer_poroelastic_small_def`; the deviation is 0.0375 $p_0$.
- `test_cryer_poroelastic_large_def_compressible`: large-deformation model, otherwise as
  `test_cryer_poroelastic_small_def_compressible`; the deviation is 0.043 $p_0$.

## Running the benchmark

In `test/poromechanics/poroelastic_large_deformations/cryer` of the build directory, `run_benchmark.py`
generates `sphere_eighth.msh` with gmsh if it does not exist, runs the cases of the table and plots them
with `plot_cryer.py` (requires numpy and matplotlib); `plot_cryer.py` writes the figure above as
`cryer_overview.png`:

```bash
python3 run_benchmark.py <build directory>
```

A single case of each model, in `test/poromechanics` of the build directory:

```bash
gmsh -3 -format msh2 sphere_eighth.geo -o sphere_eighth.msh
cd poroelastic_large_deformations/cryer
./test_cryer_poroelastic_large_def params.input -Problem.PressureLoad 100.0 -SpatialParams.PermeabilityModel Constant
cd ../../poroelastic/cryer
./test_cryer_poroelastic_small_def params.input
```

Each run writes the centre time series to `cryer_center_pressure.csv` in its directory.
