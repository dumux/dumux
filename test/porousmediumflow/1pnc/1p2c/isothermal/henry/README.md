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
velocity-dependent dispersion.

For each component $\kappa\in\{\text{solvent},\,\text{solute}\}$ we solve a mass balance,
which gives two coupled equations:

$$\frac{\partial(\phi\varrho X^\kappa)}{\partial t} +
\nabla\cdot\left(\varrho X^\kappa \mathbf{v} - \varrho D^\kappa_\text{pm}\nabla
X^\kappa\right) = q, \qquad
\mathbf{v}=-\frac{\mathbf{K}}{\mu}\left(\nabla p - \varrho\,\mathbf{g}\right)$$

where

- $\phi$ is the porosity,
- $\varrho$ is the mass density, $\varrho=\varrho_0+(\varrho_1-\varrho_0)\,
  X^\text{solute}/X_\text{sw}$, with $\varrho_0$/$\varrho_1$ the fresh-/seawater
  reference densities and $X_\text{sw}=0.035$ the seawater mass
  fraction; $c:=X^\text{solute}/X_\text{sw}\in[0,1]$ is the paper's dimensionless
  concentration ($c=0$ freshwater, $c=1$ seawater),
- $\mathbf{K}$ is the permeability tensor,
- $\mu$ is the dynamic viscosity,
- $\mathbf{v}$ is the Darcy velocity,
- $X^\kappa$ is the mass fraction of component $\kappa$
- $D^\kappa_\text{pm}$ is component $\kappa$'s diffusion/dispersion tensor,
  $D^\text{solute}_\text{pm}=\phi\,\tau D_m\mathbf{I}+\mathbf{D}$, with $D_m$ the
  molecular diffusion coefficient, $\tau$ the tortuosity. To match eq. 3 (has only $\varepsilon D_m$), we set to $\tau=1$ here so this term reduces
  to $\phi D_m$. Additionally,  $\mathbf{D}$ is Scheidegger's
  velocity-dependent dispersion tensor,
  $$\mathbf{D}=(\alpha_L-\alpha_T)\,\frac{\mathbf{v}\,\mathbf{v}^T}{|\mathbf{v}|}
  +\alpha_T\,|\mathbf{v}|\,\mathbf{I},$$
  with longitudinal/transverse dispersivities $\alpha_L$, $\alpha_T$,
- $q$ is a source/sink term, zero everywhere here.

Fahs et al. (2016) additionally assume the Boussinesq approximation: they drop
$\varrho$ from both balances above ($\nabla\cdot\mathbf{v}=0$ and $\phi\,\partial c/
\partial t+\mathbf{v}\cdot\nabla c-\nabla\cdot[(\phi D_m\mathbf{I}+\mathbf{D})\nabla
c]=0$), keeping it only in the buoyancy term of $\mathbf{v}$. This implementation
does not.

## Model assumptions {#henry-model-assumptions}

- **Incompressible fluid**: $\varrho$ depends only on $X^\text{solute}$.
- **Constant, salinity-independent viscosity**: $\mu=10^{-3}\,\mathrm{Pa\,s}$
- **Homogeneous, isotropic aquifer**: $\mathbf{K}=k\mathbf{I}$ and $\phi$ are uniform
  scalars not tensors or spatially varying fields.

## Boundary and initial conditions {#henry-boundary-and-initial-conditions}

- **Left** (inland): specified freshwater inflow (Neumann), $c=0$.
- **Right** (sea): Dirichlet, hydrostatic pressure using the seawater reference
  density, fixed concentration $c=1$
- **Top/bottom**: impermeable, no-flow.
- **Initial**: domain filled with seawater, hydrostatic.

## Test cases {#henry-test-cases}

Fahs et al. (2016) define three test cases, differing only in the dispersion
coefficients (their Table 1); all other physical parameters (their Table 2) are
identical:

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Domain length / depth | $\ell$ / $d$ | 3 / 1 | m |
| Freshwater recharge | $q_d$ | $6.6\times10^{-5}$ | m$^2$ s$^{-1}$ |
| Permeability | $k$ | $1.0204\times10^{-9}$ | m$^2$ |
| Porosity | $\phi$ | 0.35 | - |
| Freshwater density | $\rho_0$ | 1000 | kg m$^{-3}$ |
| Seawater density | $\rho_1$ | 1025 | kg m$^{-3}$ |
| Viscosity | $\mu$ | $10^{-3}$ | Pa s |

| Test case | Status | $D_m$ [m$^2$ s$^{-1}$] | $\alpha_L$ [m] | $\alpha_T$ [m] |
|-----------|--------|------------------------|----------------|----------------|
| 1 (classic, purely diffusive) | implemented | $18.86\times10^{-6}$ | 0 | 0 |
| 2 (velocity-dependent dispersion) | implemented | $9.43\times10^{-8}$ | 0.1 | 0.01 |
| 3 (velocity-dependent dispersion, narrow mixing zone) | **currently in the making** | $9.43\times10^{-8}$ | 0.001 | 0.0001 |

The molecular-diffusion part of $D^\text{solute}_\text{pm}$ is computed by DuMux's
`EffectiveDiffusivityModel`. Its default for @ref OnePNCModel, Millington-Quirk, gives
$\phi^{4/3}D_m$ in the fully saturated case, not the $\phi D_m$ of Fahs et al. (2016).
`properties.hh` therefore uses `DiffusivityConstantTortuosity` ($\phi\,\tau D_m$) with
`SpatialParams.Tortuosity = 1`, so that
$D^\text{solute}_\text{pm}=\phi D_m\mathbf{I}+\mathbf{D}$ exactly as above. With the
default, molecular diffusion would be about 30% too small
($\phi^{4/3}/\phi=0.35^{1/3}\approx0.70$).

## Coarse Fixed-grid tests (CI Test) {#henry-fixed-grid-tests}

### Setup {#henry-setup}

The domain is discretized with a structured @ref Dune::YaspGrid, 120x40 cells for
both Test Cases 1 and 2. The @ref OnePNCModel with @ref BoxDiscretization is used. Time integration uses a fixed
number of equally sized time steps to reach steady state (see `params_case1.input` / `params_case2.input`).

### Validation {#henry-validation}

Each test case is checked in two ways. The simulation runs only once per case; the
second check reuses its output.

1. **Regression check** (`test_1p2c_henry_fahs_case1_box_regression` /
   `test_1p2c_henry_fahs_case2_box_regression`): runs the simulation and compares the
   complete VTU output field by field with a stored reference result
   (`test/references/test_1p2c_henry_fahs_case<N>_box-reference.vtu`). This detects any
   unintended change of the solution. The reference is an earlier result of this code,
   so the check guards against changes, not against errors.
2. **Comparison with Fahs et al. (2016)** (`test_1p2c_henry_fahs_case1_box` /
   `test_1p2c_henry_fahs_case2_box`): `validate_fahs2016.py` takes the output of the
   regression run, extracts the $x$-positions of the 10/50/90% isochlors
   ($c=0.1,0.5,0.9$) at each depth $Z$ listed in Tables D1/D2 of the paper, and compares
   them with the table values. The test passes if the maximum relative error is below
   0.05. This tolerance is deliberately generous, since the test runs on a coarse grid
   to stay fast enough for CI.

### Results {#henry-results}

To run a test case (both the regression-check and the table-validation target) against
the corresponding table/reference:

```bash
cd <build-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry
ctest -R test_1p2c_henry_fahs_case1_box       # Test Case 1, vs. Table D1 + regression reference
ctest -R test_1p2c_henry_fahs_case2_box       # Test Case 2, vs. Table D2 + regression reference
```

@note Both targets need Python packages: the regression checks compare the output with
[`fieldcompare`](https://pypi.org/project/fieldcompare/) (through `dumux_runtest.py`),
and `validate_fahs2016.py` reads it with `meshio`. Without them the tests fail. The
easiest way is a Python virtual environment with the DuMux requirements (they include
`fieldcompare[all]`, which also brings `meshio`), set up once in the `dumux` source
directory and activated before running `ctest`:

```bash
python -m venv dumux_venv
source dumux_venv/bin/activate
pip install -r requirements.txt
```

To reproduce an animated view of the transient approach to steady state, together
with a final comparison against the literature isochlors, first run the two `ctest`
invocations above (this already runs both simulations to completion), then reuse
their output instead of re-running the simulations (from the same build directory;
requires `pyvista`, installable via `pip install pyvista` into the `dumux_venv` used
for the fuzzy/validation tooling):

```bash
python3 <source-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry/post_processing.py \
  test_1p2c_henry_fahs_case1_box.pvd test_1p2c_henry_fahs_case2_box.pvd --out henry_combined.gif
```

This produces **`henry_combined.gif`**: Test Case 1
(top) and Test Case 2 (bottom) stacked into a single animation, both sampled at the
same simulated times over $t\in[0,1]$. Each panel draws the simulated 10/50/90% isochlors as solid contour lines with the
literature Table D1/D2 points overlaid as markers. An output name ending in `.png`, e.g.
`--out henry_combined_final.png`, gives a static image of just the final-time ($t=1$ d)
isochlors against the tables instead.

![Henry problem, Test Cases 1 and 2](henry_combined.gif)

## Adaptive-grid benchmark (manual) {#henry-adaptive-benchmark}

In addition to the tests above, `main_benchmark.cc` provides two executables,
`test_1p2c_henry_case1_benchmark` and `test_1p2c_henry_case2_benchmark`, that solve the
same two test cases on an adaptive grid. They are not part of the test suite and are
run by hand. They require `dune-alugrid` and are only built on request, e.g. with
`make test_1p2c_henry_case1_benchmark`.

@note The adaptive runs take much longer than the tests above. On a single core of a
laptop (Intel i7-1260P, release build), Test Case 1 took about 2 minutes and Test
Case 2 about 50 minutes, compared to about 20 s and 1 minute on the fixed grid.

The grid starts at the same 120x40 resolution as above and is refined and coarsened
during the run so that it follows the moving saltwater/freshwater front (see
`adaptive/gridadaptindicator.hh`; the settings are in `params_benchmark_case1.input` /
`params_benchmark_case2.input`). Unlike the rectangular `YaspGrid` above, it consists of
triangles: each of the 120x40 rectangles is split in two, giving 9600 cells initially.
Triangles can be refined locally without hanging nodes, which the box scheme cannot
handle. The linear systems are solved with UMFPack, a direct solver for a single
process, so that the results reflect the adaptive grid and not the linear solver.

`post_processing.py` plots the output the same way as above. With `--grid`, it shows the
mesh instead, with each edge colored by its concentration, so you can see the
refinement follow the front:

```bash
cd <build-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry
./test_1p2c_henry_case1_benchmark params_benchmark_case1.input -Problem.Name adaptive_case1
./test_1p2c_henry_case2_benchmark params_benchmark_case2.input -Problem.Name adaptive_case2
python3 <source-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry/post_processing.py \
  adaptive_case1.pvd adaptive_case2.pvd --out henry_adaptive_solution.gif
python3 <source-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry/post_processing.py \
  adaptive_case1.pvd adaptive_case2.pvd --grid --out henry_adaptive_grid.gif
```

On the adaptive grid, the results match the benchmark much better than on the fixed
grid. The maximum relative error of the isochlor positions drops from 0.0228 to 0.0107
for Test Case 1 and from 0.0417 to 0.0129 for Test Case 2 (same check as
`validate_fahs2016.py` above). This is also visible in the plots: in the top right of
the domain, the simulated $c=0.1$ isochlor of Test Case 2 now runs through the
reference points instead of next to them.

![Henry problem, adaptive refinement, solution fit](henry_adaptive_solution.gif)

The mesh plot shows whether the adaptivity works as intended: the grid is refined
where the concentration changes quickly, along the mixing zone between fresh and salt
water, and stays at the coarse starting resolution away from the front, most visibly in
the freshwater region on the left. As the front moves, the fine region moves with it.

![Henry problem, adaptive refinement, mesh colored by concentration](henry_adaptive_grid.gif)
