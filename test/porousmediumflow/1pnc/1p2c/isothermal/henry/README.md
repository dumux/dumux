# Benchmark: Henry Saltwater Intrusion Problem {#benchmark-henry}

**Problem Description**

The Henry problem (Henry 1964 @cite Henry1964) is a classic benchmark for
density-driven groundwater flow and solute transport, describing seawater intrusion
into a confined, homogeneous coastal aquifer. Freshwater is injected at a fixed rate
along the left (inland) boundary; the right (seaward) boundary is in hydrostatic
contact with seawater at fixed salinity. The resulting steady state is a saltwater
wedge intruding along the bottom of the aquifer beneath an outflowing freshwater lens.

Considered is the benchmark of Fahs et al. (2016) @cite Fahs2016, itself a
re-derivation of Henry (1964) @cite Henry1964 with a much larger number of Fourier
modes and an extension to velocity-dependent (Scheidegger) dispersion.

Implemented is one equation, solved for each component $\kappa\in\{\text{solvent},\,
\text{solute}\}$ (two coupled PDEs in total), the component mass balance

$$\frac{\partial(\phi\varrho X^\kappa)}{\partial t} -
\nabla\cdot\left(\varrho X^\kappa \mathbf{v} + \varrho D^\kappa_\text{pm}\nabla
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

**Model Assumptions**

- **Incompressible fluid**: $\varrho$ depends only on $X^\text{solute}$.
- **Constant, salinity-independent viscosity**: $\mu=10^{-3}\,\mathrm{Pa\,s}$
- **Homogeneous, isotropic aquifer**: $\mathbf{K}=k\mathbf{I}$ and $\phi$ are uniform
  scalars not tensors or spatially varying fields.

**Boundary and Initial Conditions**

- **Left** (inland): specified freshwater inflow (Neumann), $c=0$.
- **Right** (sea): Dirichlet, hydrostatic pressure using the seawater reference
  density, fixed concentration $c=1$
- **Top/bottom**: impermeable, no-flow.
- **Initial**: domain filled with seawater, hydrostatic.

**Test Cases**

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

Test Case 1 is purely diffusive: $\mathbf{D}=\mathbf{0}$, i.e. no dispersion tensor at
all (dispersion is disabled entirely). Test Case 2 has $\mathbf{D}$ nonzero, enabling
the Scheidegger dispersion tensor already built into DuMux
(`Dumux::ScheideggersDispersionTensor`,`dumux/material/fluidmatrixinteractions/dispersiontensors/scheidegger.hh`)
via `EnableCompositionalDispersion` and `spatialparams.hh`'s `dispersionAlphas()`.
Test Case 3 uses the same machinery with smaller dispersivities.

**Setup**

The domain is discretized with a structured @ref Dune::YaspGrid, 240x80 cells for
both Test Cases 1 and 2 (the paper uses a finer mesh for Test Case 2, not replicated
here -- see `params_case2.input`). The @ref
OnePNCModel with @ref BoxDiscretization is used. Time integration uses a fixed
number of equally sized time steps, run well past the time the paper
reports for each case to reach steady state (see `params.input` / `params_case2.input`
for the exact, and explicitly *not* empirically verified, margins chosen).

DuMux's default `EffectiveDiffusivityModel` for @ref OnePNCModel is Millington-Quirk
($D_\mathrm{eff}=D_m\phi^{1/3}$), which does not match the paper's transport
equation above ($\phi D_m$, linear in porosity, no separate tortuosity
reduction). `properties.hh` overrides this to `DiffusivityConstantTortuosity` with
`SpatialParams.Tortuosity = 1`, giving $D_\mathrm{eff}=\phi D_m$ exactly. Left
at the DuMux default ($\phi=0.35$), this would apply roughly $2\times$ too much
molecular diffusion ($0.35^{1/3}/0.35\approx2.01$).

**Validation**

Fahs et al. (2016) digitized their converged semianalytical isochlor
positions (Appendix D, Tables D1/D2/D3). The test extracts the simulated 10/50/90% isochlor
($c=0.1,0.5,0.9$) $x$-positions at each tabulated depth $Z$ by linear interpolation
of the concentration field, and compares them directly against those table values.

**Results**

To run a test case and validate it against the corresponding table:

```bash
cd <build-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry
ctest -R test_1p2c_henry_fahs_box       # Test Case 1, vs. Table D1
ctest -R test_1p2c_henry_fahs_case2_box # Test Case 2, vs. Table D2
```

To reproduce an animated view of the transient approach to steady state, together
with a final comparison against the literature isochlors, first run the two `ctest`
invocations above (this already runs both simulations to completion), then reuse
their output instead of re-running the simulations (from the same build directory;
requires `pyvista`, installable via `pip install pyvista` into the `dumux_venv` used
for the fuzzy/validation tooling):

```bash
python3 <source-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry/make_gif.py --combined --skip-run
```

This produces **`henry_combined.gif`**: Test Case 1
(top) and Test Case 2 (bottom) stacked into a single animation, both sampled at the
same simulated times over $t\in[0,1]$ d, i.e. the full run (`TimeLoop.TEnd` is 1 d for
both cases; both are converged well before that). Rendering both cases into one image is
deliberate: two separate GIF files would each start playing on their own
load/decode schedule and drift out of sync in a browser regardless of matching
frame timing, whereas a single combined image is in sync by construction. Each
panel draws the simulated 10/50/90% isochlors as solid contour lines with the
literature Table D1/D2 points overlaid as markers. A static
**`henry_combined_final.png`** comparing the final-time ($t=1$ d) isochlors against
the tables is saved the same way.

`--case 1`/`--case 2` remain available for quick single-case iteration, producing
their own separate `henry_case<N>.gif`/`henry_case<N>_final.png`.

![Henry problem, Test Cases 1 and 2](henry_combined.gif)
