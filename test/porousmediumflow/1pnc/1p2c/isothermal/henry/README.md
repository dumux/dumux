# Benchmark: Henry Saltwater Intrusion Problem {#benchmark-henry}

The problem of Henry (1964) @cite Henry1964 describes seawater intrusion into a confined coastal aquifer:
freshwater flows towards the sea, and the denser seawater forms a wedge along the bottom of the aquifer.
Fahs et al. (2016) @cite Fahs2016 computed semianalytical solutions (Fourier-Galerkin) of the steady state for
molecular diffusion only (their test case 1) and with velocity-dependent dispersion (test case 2) and tabulated the
positions of the 10, 50 and 90 % isochlors (their Tables D1 and D2). The tests compare the steady state computed with
DuMux with these tables.

## Problem description {#henry-problem-description}

The aquifer is a rectangle of length $\ell = 3\,\mathrm{m}$ and depth $d = 1\,\mathrm{m}$. With the Boussinesq
approximation, the Darcy velocity $\mathbf{q}$ and the concentration $c$ (1 for seawater) satisfy (eqs. 1-5 of
Fahs et al. (2016) @cite Fahs2016)

$$
\nabla \cdot \mathbf{q} = 0, \qquad
\mathbf{q} = -\frac{k}{\mu} \left( \nabla p - (\varrho_0 + \Delta\varrho\, c)\, \mathbf{g} \right),
$$

$$
\phi \frac{\partial c}{\partial t} + \mathbf{q} \cdot \nabla c - \nabla \cdot \left( (\phi D_m \mathbf{I} + \mathbf{D}) \nabla c \right) = 0,
\qquad
\mathbf{D} = (\alpha_L - \alpha_T) \frac{\mathbf{q} \mathbf{q}^T}{|\mathbf{q}|} + \alpha_T |\mathbf{q}| \mathbf{I},
$$

with the permeability $k$, the viscosity $\mu$, the freshwater density $\varrho_0$, the density difference
$\Delta\varrho$ between seawater and freshwater, the porosity $\phi$, the molecular diffusion coefficient $D_m$ and the
longitudinal and transverse dispersivities $\alpha_L$ and $\alpha_T$.

Boundary and initial conditions:

- Sea ($x = \ell$): $c = 1$ and hydrostatic seawater pressure, $p = p_0 + (\varrho_0 + \Delta\varrho) g (d - z)$.
- Inland ($x = 0$): $c = 0$ and no vertical velocity, i.e. hydrostatic freshwater pressure,
  $p = p_\mathrm{in} + \varrho_0 g (d - z)$. Since $\nabla \cdot \mathbf{q} = 0$ and top and bottom are impermeable,
  the inflow per unit width is the same through every vertical cross section. Integrating $q_x$ over the domain gives
  $q_d = \frac{k d}{\mu \ell} \left( p_\mathrm{in} - p_0 - \frac{1}{2} \Delta\varrho\, g d \right)$ for any
  concentration field, so the prescribed inflow $q_d$ is obtained with
  $p_\mathrm{in} = p_0 + \frac{q_d \mu \ell}{k d} + \frac{1}{2} \Delta\varrho\, g d$.
- Top and bottom: no flow.
- Initially, the aquifer is filled with seawater at rest. The simulation runs for one day, after which the
  solution is steady.

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Freshwater inflow per unit width | $q_d$ | $6.6 \times 10^{-5}\,\mathrm{m^2/s}$ |
| Permeability | $k$ | $1.0204 \times 10^{-9}\,\mathrm{m^2}$ |
| Porosity | $\phi$ | $0.35$ |
| Freshwater density | $\varrho_0$ | $1000\,\mathrm{kg/m^3}$ |
| Density difference | $\Delta\varrho$ | $25\,\mathrm{kg/m^3}$ |
| Viscosity | $\mu$ | $10^{-3}\,\mathrm{Pa\,s}$ |
| Gravitational acceleration | $g$ | $9.81\,\mathrm{m/s^2}$ |

| Case | $D_m$ | $\alpha_L$ | $\alpha_T$ |
|------|-------|------------|------------|
| Molecular diffusion | $1.886 \times 10^{-5}\,\mathrm{m^2/s}$ | 0 | 0 |
| Velocity-dependent dispersion | $9.43 \times 10^{-8}\,\mathrm{m^2/s}$ | $0.1\,\mathrm{m}$ | $0.01\,\mathrm{m}$ |

## Model and discretization {#henry-model}

The tests use the @ref OnePNCModel with mass fractions, where the solute mass fraction is $c$. Two components with
equal molar masses give a constant fluid density $\varrho_0$, and the Darcy flux `CVFEBoussinesqDarcysLaw` adds the
buoyancy $\Delta\varrho\, c$ in the gravity term (`problem.buoyantDensity`). With a constant density, the mass
balances of the two components reduce to the equations above. `DiffusivityConstantTortuosity` with tortuosity 1
gives $\phi D_m$, and `ScheideggersDispersionTensor` gives $\mathbf{D}$.

The grid has $120 \times 40$ cells and is discretized with the @ref BoxDiscretization. The advective fluxes use the
arithmetic mean of the upstream and downstream values (`Flux.UpwindWeight = 0.5`). The cell Péclet number
$h |\mathbf{q}| / (\phi D_m + \alpha_L |\mathbf{q}|)$ is at most 2.0 with molecular diffusion only (at the outflow)
and 0.25 with dispersion, and the concentration stays within $[0, 1]$. The implicit Euler method runs with time
steps of at most one hour.

## Results {#henry-results}

The top row of the figure shows the steady concentration field with the computed 10, 50 and 90 % isochlors (lines),
all isochlor positions of Fahs et al. (2016) @cite Fahs2016 (circles) and the midline $z = 0.5\,\mathrm{m}$
(dashed), along which the concentration is plotted over time below.

![Steady concentration field with isochlors and reference positions (top) and the concentration along the midline over time (light to dark).](henry_profiles.png)

Maximum (and root mean square) distance between the computed isochlors and the tabulated positions, 57 points
per case:

| Grid | Molecular diffusion | Velocity-dependent dispersion |
|------|-------------|-------------|
| $120 \times 40$ | 5.2 mm (1.0 mm) | 7.0 mm (2.6 mm) |
| $240 \times 80$ | 5.9 mm (1.1 mm) | 4.9 mm (1.4 mm) |
| $480 \times 160$ | 6.6 mm (1.1 mm) | 5.4 mm (1.6 mm) |

Positions $x$ (in m) of the isochlors on the midline $z = 0.5\,\mathrm{m}$:

| Case | Isochlor | Fahs et al. (2016) @cite Fahs2016 | $120 \times 40$ | $240 \times 80$ | $480 \times 160$ |
|------|----------|--------------------|------------------|------------------|-------------------|
| Molecular diffusion | 10 % | 2.291 | 2.291 | 2.291 | 2.291 |
|  | 50 % | 2.748 | 2.748 | 2.748 | 2.748 |
|  | 90 % | 2.964 | 2.963 | 2.963 | 2.963 |
| Velocity-dependent dispersion | 10 % | 2.438 | 2.440 | 2.439 | 2.438 |
|  | 50 % | 2.714 | 2.712 | 2.713 | 2.713 |
|  | 90 % | 2.881 | 2.883 | 2.881 | 2.880 |

## Running the benchmark {#henry-running}

The tests `test_1p2c_henry_diffusion_box` and `test_1p2c_henry_dispersion_box` fail if the maximum distance exceeds
$1\,\mathrm{cm}$. The figure is created from their output:

```bash
cd <build-dir>/test/porousmediumflow/1pnc/1p2c/isothermal/henry
make test_1p2c_henry_diffusion_box
ctest -R test_1p2c_henry_
python3 plot_profiles.py test_1p2c_henry_diffusion_box.pvd test_1p2c_henry_dispersion_box.pvd
```

Test case 3 of Fahs et al. (2016) @cite Fahs2016 ($\alpha_L = 1\,\mathrm{mm}$) needs cells smaller than about
$2\,\mathrm{mm}$ for a cell Péclet number below 2 and is not part of the tests.
