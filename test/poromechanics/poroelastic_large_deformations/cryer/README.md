# Cryer's Sphere Consolidation {#benchmark-cryer}

A fluid-saturated porous sphere is loaded on its drained surface. The pore pressure at the centre
first rises above the load and then decays (Mandel–Cryer effect). Cryer @cite Cryer1963 solved the
problem in linear poroelasticity. The test solves it with the large-deformation poroelastic model
(@ref PoroElasticLargeDef): with a small load the centre pressure is compared with Cryer's solution,
with larger loads the test shows the effect of finite deformation and of a permeability that depends on
the deformation.

## Problem

The sphere has the radius \f$ R_0 = 1 \f$ m; one octant is modelled.
- Sphere surface: the dead pressure load \f$ p_0 \f$ (nominal traction \f$ -p_0 \mathbf{N} \f$ on the reference
  surface) and the drained condition \f$ p_f = 0 \f$.
- Symmetry planes \f$ x = 0 \f$, \f$ y = 0 \f$, \f$ z = 0 \f$: zero normal displacement and no flow.

![Section z = 0 of the simulated octant: load and drainage on the sphere surface, rollers on the symmetry planes x = 0 and y = 0](cryer_setup.svg){html: width=50%}

The skeleton is a Neo-Hookean material with volumetric–isochoric split, isochoric stress
\f$ \mathbf{P}_\text{iso} = G J^{-2/3}(\mathbf{F} - \tfrac{1}{3}\operatorname{tr}(\mathbf{C})\, \mathbf{F}^{-\top}) \f$ and
solid pressure \f$ p_s = \frac{K}{2}(J - J^{-1}) \f$. The model adds \f$ (p_s - \alpha_B p_f) J \mathbf{F}^{-\top} \f$ to the
stress and has the storage \f$ \alpha_B \partial_t J + S_p \partial_t p_f \f$; see @ref PoroElasticLargeDef for the
equations.

| Parameter | Symbol | Value |
| :--- | :--- | ---: |
| shear modulus | \f$ G \f$ | 1.5 MPa |
| drained bulk modulus | \f$ K \f$ | 1.0 MPa |
| Poisson's ratio of the linearized skeleton | \f$ \nu \f$ | 0 |
| Biot coefficient | \f$ \alpha_B \f$ | 1 |
| storage coefficient | \f$ S_p \f$ | 0 |
| initial porosity | \f$ \phi_0 \f$ | 0.5 |
| initial permeability | \f$ \kappa_0 \f$ | \f$ 10^{-12} \f$ m² |
| fluid viscosity | \f$ \mu_f \f$ | \f$ 10^{-3} \f$ Pa s |
| surface load | \f$ p_0 \f$ | 100 Pa, 0.25 MPa or 0.5 MPa (\f$ p_0/K = 10^{-4} \f$, 0.25, 0.5) |

The permeability is either constant, \f$ \kappa = \kappa_0 \f$, or follows the Kozeny–Carman law
\f$ \kappa(J) = \frac{d_p^2}{180} \frac{\phi^3}{(1 - \phi)^2} \f$ with the porosity \f$ \phi(J) = 1 - (1 - \phi_0)/J \f$ and the
particle diameter \f$ d_p \f$ chosen such that \f$ \kappa(1) = \kappa_0 \f$ (`SpatialParams.PermeabilityModel`).

## Cryer's solution

In linear poroelasticity, the pore pressure is
\f[
\frac{p(r,t)}{p_0} = 2mB \sum_{i=1}^{\infty}
   \frac{\sin z_i - z_i \operatorname{sinc}(z_i r/R_0)}{m z_i \cos z_i + (2m-1)\sin z_i}
   \exp\left(-\frac{z_i^2 c_v t}{\beta R_0^2}\right), \qquad \operatorname{sinc}(x) = \frac{\sin x}{x},
\f]
where the \f$ z_i \f$ are the positive roots of \f$ (1 - m z^2)\sin z - z\cos z = 0 \f$ and
\f[
m = \frac{K + \tfrac{4}{3}G}{4G}\left(1 + \frac{K S_p}{\alpha_B^2}\right), \qquad
c_v = \frac{\kappa}{\mu_f}\left(K + \tfrac{4}{3}G\right), \qquad
\beta = \alpha_B^2 + \left(K + \tfrac{4}{3}G\right) S_p, \qquad
B = \frac{\alpha_B}{\alpha_B^2 + K S_p}.
\f]
\f$ B \f$ is Skempton's coefficient, the uniform pore pressure \f$ B p_0 \f$ right after loading. With the
parameters above, \f$ m = 0.5 \f$, \f$ c_v = 3 \times 10^{-3} \f$ m²/s, \f$ \beta = 1 \f$, \f$ B = 1 \f$, and the centre
pressure has its maximum 1.575 \f$ p_0 \f$ at \f$ t c_v/R_0^2 = 0.0625 \f$. `test/poromechanics/cryer_analytical.hh` implements the series.

## Discretization

The displacement is quadratic with the hybrid control-volume / finite-element scheme
(@ref PQ2Discretization), solid and fluid pressure are linear with the box scheme
(@ref BoxDiscretization), coupled as three subdomains (@ref MultiDomain). The tests use the mesh
`test/poromechanics/sphere_eighth_coarse.msh` (223 vertices, 670 tetrahedra). The benchmark runs use the
mesh `sphere_eighth.msh` that gmsh generates from `test/poromechanics/sphere_eighth.geo` with the maximum
element size
0.095 m; with gmsh 4.12.1 it has 854 vertices and 3238 tetrahedra. Time integration uses the multi-stage time stepper with the
third-order diagonally implicit Runge–Kutta method DIRK3 (`TimeLoop.Scheme`, implicit Euler is also
available). The time step starts at 1 s and grows by the factor 1.3 per step up to 20 s, which reaches
\f$ R_0^2/c_v = 333.33 \f$ s in 25 steps.

## Results

Centre pressure \f$ p/p_0 \f$ and volume ratio \f$ J \f$ at the centre, mesh `sphere_eighth.msh` (gmsh 4.12.1), DIRK3,
25 steps:

| permeability | \f$ p_0/K \f$ | maximum of \f$ p/p_0 \f$ (at \f$ t c_v/R_0^2 \f$) | \f$ J \f$ at \f$ t c_v/R_0^2 = 1 \f$ | \f$ p/p_0 \f$ at \f$ t c_v/R_0^2 = 1 \f$ |
| :--- | ---: | ---: | ---: | ---: |
| constant | \f$ 10^{-4} \f$ | 1.583 (0.072) | 0.9999 | 0.031 |
| constant | 0.25 | 1.710 (0.072) | 0.743 | 0.028 |
| constant | 0.5 | 1.865 (0.072) | 0.473 | 0.017 |
| Kozeny–Carman | \f$ 10^{-4} \f$ | 1.583 (0.072) | 0.9999 | 0.031 |
| Kozeny–Carman | 0.25 | 1.587 (0.072) | 0.802 | 0.506 |
| Kozeny–Carman | 0.5 | 1.566 (0.096) | 0.738 | 1.059 |

With the small load, the centre pressure deviates from Cryer's solution by at most 0.018 \f$ p_0 \f$ at the
25 time steps; at \f$ t c_v/R_0^2 = 0.072 \f$, Cryer's solution is 1.568 \f$ p_0 \f$. With the larger loads, the
overshoot grows for a constant permeability. With the Kozeny–Carman law, the compaction reduces the
permeability, for \f$ J = 0.738 \f$ to 15 % of \f$ \kappa_0 \f$, and the centre has not drained at
\f$ t c_v/R_0^2 = 1 \f$.

![Centre pressure (top) and volume ratio at the centre (bottom) for the three loads, constant (left) and Kozeny–Carman permeability (right); dashed: Cryer's solution](cryer_overview.png)

## Tests

- `test_cryer_poroelastic_large_def`: coarse mesh, implicit Euler, five steps of 1 s, small load. The
  test fails if the centre pressure deviates from Cryer's solution by more than 0.07 \f$ p_0 \f$
  (`Problem.CenterPressureTolerance`); the deviation is 0.0375 \f$ p_0 \f$.
- `test_cryer_poroelastic_large_def_compressible`: as above with \f$ \alpha_B = 0.8 \f$ and
  \f$ S_p = 5 \times 10^{-7} \f$ 1/Pa (\f$ B = 0.70 \f$), to 180 s with growing time steps; the deviation is
  0.043 \f$ p_0 \f$.

## Running the benchmark

In `test/poromechanics/poroelastic_large_deformations/cryer` of the build directory, `run_benchmark.py`
generates `sphere_eighth.msh` with gmsh if it does not exist, runs the six cases of the table and plots
them with `plot_cryer.py` (requires numpy and matplotlib); `plot_cryer.py` writes the figure above as
`cryer_overview.png`:

```bash
python3 run_benchmark.py <build directory>
```

A single case, in `test/poromechanics` of the build directory:

```bash
gmsh -3 -format msh2 sphere_eighth.geo -o sphere_eighth.msh
cd poroelastic_large_deformations/cryer
./test_cryer_poroelastic_large_def params.input -Problem.PressureLoad 100.0 -SpatialParams.PermeabilityModel Constant
python3 plot_cryer.py
```

Each run writes the centre time series to `cryer_center_pressure.csv`.
