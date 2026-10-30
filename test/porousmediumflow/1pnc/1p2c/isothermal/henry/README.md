# Benchmark: Henry Saltwater Intrusion Problem {#benchmark-henry}

## Problem description {#henry-problem-description}

The Henry problem (Henry 1964 @cite Henry1964) is a classic benchmark for
density-driven groundwater flow and solute transport, describing seawater intrusion
into a confined, homogeneous coastal aquifer. Freshwater is injected at a fixed rate
along the left (inland) boundary; the right (seaward) boundary is in hydrostatic
contact with seawater at fixed salinity. The resulting steady state is a saltwater
wedge intruding along the bottom of the aquifer beneath an outflowing freshwater lens.

We consider the benchmark of Fahs et al. (2016) @cite Fahs2016, a re-derivation of
Henry (1964) @cite Henry1964 with higher accuracy and an extension to
velocity-dependent dispersion. They tabulate the positions of the 10, 50 and 90 %
isochlors of the steady state (their Tables D1 and D2), which the tests compare with.

Fahs et al. (2016) use the Boussinesq approximation: the density $\varrho$ is dropped from
both balances and kept only in the buoyancy term of the Darcy velocity. With the
Darcy velocity $\mathbf{v}$ and the dimensionless concentration $c\in[0,1]$
($c=0$ freshwater, $c=1$ seawater), the problem reads (eqs. 1-5 of Fahs et al. (2016) @cite Fahs2016)

$$\nabla\cdot\mathbf{v}=0, \qquad
\mathbf{v}=-\frac{\mathbf{K}}{\mu}\left(\nabla p - (\varrho_0+\Delta\varrho\,c)\,\mathbf{g}\right),$$

$$\phi\,\frac{\partial c}{\partial t}+\mathbf{v}\cdot\nabla c-\nabla\cdot\left[(\phi D_m\mathbf{I}+\mathbf{D})\nabla c\right]=0,
\qquad
\mathbf{D}=(\alpha_L-\alpha_T)\,\frac{\mathbf{v}\,\mathbf{v}^T}{|\mathbf{v}|}+\alpha_T\,|\mathbf{v}|\,\mathbf{I},$$

where

- $\phi$ is the porosity,
- $\varrho_0$ is the freshwater density and $\Delta\varrho$ the density difference between
  seawater and freshwater,
- $\mathbf{K}=k\mathbf{I}$ is the permeability tensor,
- $\mu$ is the dynamic viscosity,
- $D_m$ is the molecular diffusion coefficient,
- $\mathbf{D}$ is Scheidegger's velocity-dependent dispersion tensor with the
  longitudinal/transverse dispersivities $\alpha_L$, $\alpha_T$.

## Model assumptions {#henry-model-assumptions}

- **Boussinesq approximation**: the fluid density is constant, only the buoyancy
  term depends on the concentration.
- **Constant, salinity-independent viscosity**: $\mu=10^{-3}\,\mathrm{Pa\,s}$
- **Homogeneous, isotropic aquifer**: $\mathbf{K}=k\mathbf{I}$ and $\phi$ are uniform
  scalars, not tensors or spatially varying fields.

## Boundary and initial conditions {#henry-boundary-and-initial-conditions}

- **Right** (sea): Dirichlet, $c=1$ and hydrostatic seawater pressure
  $p=p_0+(\varrho_0+\Delta\varrho)\,g\,(d-z)$.
- **Left** (inland): uniform freshwater inflow $q_d/d$ (Neumann), with zero total solute flux. The solute
  concentration is not prescribed: in the steady state the saltwater wedge is far from the left boundary, and
  $c$ is below $10^{-9}$ there in the simulations.
- **Top/bottom**: impermeable, no-flow.
- **Initial**: domain filled with seawater at rest.

## Test cases {#henry-test-cases}

Fahs et al. (2016) define three test cases, differing only in the diffusion and dispersion
coefficients (their Table 1); all other physical parameters (their Table 2) are identical:

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Domain length / depth | $\ell$ / $d$ | 3 / 1 | m |
| Freshwater recharge | $q_d$ | $6.6\times10^{-5}$ | m$^2$ s$^{-1}$ |
| Permeability | $k$ | $1.0204\times10^{-9}$ | m$^2$ |
| Porosity | $\phi$ | 0.35 | - |
| Freshwater density | $\varrho_0$ | 1000 | kg m$^{-3}$ |
| Density difference | $\Delta\varrho$ | 25 | kg m$^{-3}$ |
| Viscosity | $\mu$ | $10^{-3}$ | Pa s |

| Test case | Status | Input file | $D_m$<br>[m$^2$ s$^{-1}$] | $\alpha_L$<br>[m] | $\alpha_T$<br>[m] |
|-----------|--------|------------|------------------------|----------------|----------------|
| 1 (classic, purely diffusive) | implemented | `params_diffusion.input` | $18.86\times10^{-6}$ | 0 | 0 |
| 2 (velocity-dependent dispersion) | implemented | `params_dispersion.input` | $9.43\times10^{-8}$ | 0.1 | 0.01 |
| 3 (velocity-dependent dispersion, narrow mixing zone) | not implemented | - | $9.43\times10^{-8}$ | 0.001 | 0.0001 |

Test case 3 needs cells smaller than about $2\,\mathrm{mm}$ for a cell Péclet number below 2
and is not part of the tests.

## Model and discretization {#henry-model}

The tests use the @ref OnePNCModel with mass fractions, where the solute mass fraction is $c$.
Two components with equal molar masses give a constant fluid density $\varrho_0$. The Darcy flux
`CVFEBoussinesqDarcysLaw` adds the buoyancy $\Delta\varrho\,c$ in the gravity term, through the
problem function `buoyantDensity`. With a constant density, the mass balances of the two
components reduce to the equations above.

### Diffusion and dispersion {#henry-dispersion}

The solute is transported by diffusion and dispersion with the gradient $\nabla c$ and the tensor
$D_\text{eff}\,\mathbf{I}+\mathbf{D}$, where $D_\text{eff}$ is the effective diffusion coefficient. To match eq. 3 of
Fahs et al. (2016), it must be

$$D_\text{eff}=\phi\,D_m.$$

We use the model `DiffusivityConstantTortuosity`, which computes $D_\text{eff}=\phi\,S\,\tau\,D_m$. With the saturation
$S=1$ and the tortuosity $\tau=1$ (`SpatialParams.Tortuosity = 1`), this is $D_\text{eff}=\phi D_m$, the same as in
the paper. The model is set in `properties.hh`:

```cpp
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::HenryFahs>
{ using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };
```

The dispersion tensor $\mathbf{D}$ is Scheidegger's tensor of the Darcy velocity $\mathbf{v}$ (`ScheideggersDispersionTensor`,
the default dispersion model),
$$\mathbf{D}=(\alpha_L-\alpha_T)\,\frac{\mathbf{v}\,\mathbf{v}^T}{|\mathbf{v}|}+\alpha_T\,|\mathbf{v}|\,\mathbf{I}.$$
It is switched on by the property `EnableCompositionalDispersion`, and the dispersivities are read from
`SpatialParams.LongitudinalDispersivity` and `SpatialParams.TransverseDispersivity` (both 0 in `params_diffusion.input`).

The domain is discretized with a structured @ref Dune::YaspGrid, $120\times40$ cells, and the
@ref BoxDiscretization. The advective fluxes use central differences, i.e. the arithmetic
mean of the upstream and downstream values (`Flux.UpwindWeight = 0.5`). The cell Péclet number
$h|\mathbf{v}|/(\phi D_m+\alpha_L|\mathbf{v}|)$ is at most 2.0 with molecular diffusion only (at the
outflow) and 0.25 with dispersion. The implicit Euler method runs with time steps of at most one
hour, for one simulated day, after which the solution is steady. The linear systems are solved
with UMFPack.

## Validation {#henry-validation}

At the end of the simulation, `main.cc` determines the $x$-positions of the 10/50/90 % isochlors
($c=0.1,0.5,0.9$) at each height $z$ listed in Tables D1/D2 of the paper by linear interpolation
along the grid line, and prints their difference to the table. The error is thus the horizontal offset
$x_\mathrm{computed}(z)-x_\mathrm{table}(z)$ in metres at a given height $z$, not a concentration difference or the
distance to the nearest table point. The test fails if its maximum absolute value
exceeds `Problem.IsochlorTolerance`, here $1\,\mathrm{cm}$.

Maximum horizontal offset between the computed isochlors and the tabulated positions
(19 heights, 3 isochlors each, $120\times40$ cells, central differences):

| Test case | Maximum horizontal isochlor offset |
|-----------|------------------------|
| 1, molecular diffusion | 5.2 mm |
| 2, velocity-dependent dispersion | 7.0 mm |

The cell size of the CI grid is $h=25\,\mathrm{mm}$, so both maxima are a fraction of one cell.

## Convergence {#henry-convergence}

The root mean square error (RMSE) of the horizontal isochlor offsets (same 57 values per test case) is computed on grids
with 9x3 to 240x80 cells, all other settings unchanged (Neumann inflow on the left). The offset at heights between grid lines is computed from the
concentration interpolated linearly in $z$.

| Grid | $h$ [m] | Test case 1: RMSE (max) | Test case 2: RMSE (max) |
|------|---------|--------------------------|--------------------------|
| $9\times3$ | 0.333 | 165.8 mm (729 mm) | 112.1 mm (509 mm) |
| $15\times5$ | 0.200 | 79.0 mm (403 mm) | 42.1 mm (143 mm) |
| $30\times10$ | 0.100 | 23.1 mm (152 mm) | 22.1 mm (107 mm) |
| $60\times20$ | 0.050 | 2.8 mm (10.3 mm) | 8.4 mm (33 mm) |
| $120\times40$ | 0.025 | 1.05 mm (5.2 mm) | 2.6 mm (7.0 mm) |
| $240\times80$ | 0.0125 | 1.08 mm (5.9 mm) | 1.4 mm (4.9 mm) |

The study is not part of the tests. After building the test executable (`make test_1p2c_henry_diffusion_box`), run
`python3 convergence.py` in the build directory of this test. The script only runs the executable on the grids
and writes `henry_convergence.png`. It needs the Python packages of the animation below and `numpy`.

![Convergence of the isochlor positions towards the semianalytical solution](henry_convergence.png)

The coarsest grids are not in the asymptotic range: the cell Péclet number is above 2 and the saltwater wedge is
not resolved. From $60\times20$ to $120\times40$ the RMSE decreases with an observed order of 1.4 (test case 1) and
1.7 (test case 2). For test case 2 it continues to decrease on the finest grid with an order of 0.9. For test case 1
it no longer decreases, at about 1 mm. Likely causes are errors that do not depend on the grid size, such as the
time step (at most one hour) and the accuracy of the tabulated reference values; this was not tested.

## Running the benchmark {#henry-results}

The two tests are run by

```bash
cd <build-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry
make test_1p2c_henry_diffusion_box
ctest -R test_1p2c_henry_
```

They use the same executable with different input files. The output files
`test_1p2c_henry_diffusion_box.pvd` and `test_1p2c_henry_dispersion_box.pvd` are used
for the animation below. Calling `./test_1p2c_henry_diffusion_box params_dispersion.input` directly
also prints the table of the errors per height.

@note The animation needs the Python packages `pyvista`, `matplotlib` and `pillow`.
The easiest way is a Python virtual environment, set up once in the `dumux` source
directory and activated before running the script:

```bash
python -m venv dumux_venv
source dumux_venv/bin/activate
pip install pyvista matplotlib pillow
```

`post_processing.py` plots the output of both tests as an animation of the transient approach to the steady state,
sampled at equal times over $t\in[0,1]$ d. Test case 1 (top, molecular diffusion only) and test case 2 (bottom,
with dispersion) are stacked. Each panel shows the concentration, the simulated 10/50/90 % isochlors as lines,
and the literature points of Tables D1/D2 as markers. In the build directory, after running both tests:

```bash
python3 post_processing.py
```

It reads `test_1p2c_henry_diffusion_box.pvd` and `test_1p2c_henry_dispersion_box.pvd` and writes `henry_combined.gif`.

![Henry problem, test cases 1 and 2, transient approach to the steady state at 1 d, compared with Tables D1 and D2 of Fahs et al.](henry_combined.gif)
