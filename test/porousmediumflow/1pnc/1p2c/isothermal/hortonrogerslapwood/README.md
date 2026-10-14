# Onset of Convection in a Porous Layer (Horton-Rogers-Lapwood) {#benchmark-horton-rogers-lapwood}

## Growth rate of a convection cell with the Boussinesq approximation

**Problem description**

A porous layer of height $H$ between impermeable walls contains a fluid with solute mass fraction $X = 1$
at the top and $X = 0$ at the bottom. The density increases with $X$, so the heavier fluid lies on top.
The linear profile $X = z/H$ at rest is a steady state. Above a critical Rayleigh number, perturbations of
this state grow into convection cells @cite HortonRogers1945 @cite Lapwood1948. The benchmark compares the
growth rate of a single convection cell with linear stability theory. It tests Darcy's law with the
Boussinesq approximation, `CVFEBoussinesqDarcysLaw`.

With the Boussinesq approximation, the density varies only in the gravity term of Darcy's law,

$$\nabla \cdot \mathbf{v} = 0,
\qquad
\mathbf{v} = -\frac{K}{\mu} \left( \nabla p - (\varrho + \Delta\varrho X) \mathbf{g} \right),
\qquad
\phi \frac{\partial X}{\partial t} + \nabla \cdot \left( X \mathbf{v} - D_\mathrm{eff} \nabla X \right) = 0,$$

where $\mathbf{v}$ is the Darcy velocity, $p$ the pressure, $K$ the permeability, $\mu$ the viscosity,
$\varrho$ the constant reference density of the fluid, $\Delta\varrho$ the density increase at $X = 1$,
$\mathbf{g}$ the gravitational acceleration with magnitude $g$, $\phi$ the porosity and
$D_\mathrm{eff} = \phi \tau D$ the effective diffusion coefficient with tortuosity $\tau$ and molecular
diffusion coefficient $D$. The problem returns the density deviation $\Delta\varrho X$ with
`buoyantDensity(volVars, phaseIdx)`; the fluid system returns the reference density $\varrho$, which is also used
in the storage and advection terms.

**Linear stability theory**

The only parameter of the problem is the Rayleigh number

$$Ra = \frac{K \Delta\varrho g H}{\mu D_\mathrm{eff}}.$$

With lengths in units of $H$ and time in units of $\phi H^2 / D_\mathrm{eff}$, a perturbation
$X - z = A(t) \cos(k x) \sin(\pi z)$ of the linear profile with wavenumber $k$ satisfies the linearized
equations with $A \propto e^{\sigma t}$ and the growth rate

$$\sigma = Ra \frac{k^2}{k^2 + \pi^2} - \left( k^2 + \pi^2 \right).$$

Perturbations grow for $Ra > (k^2 + \pi^2)^2 / k^2$. The smallest such Rayleigh number is
$Ra_c = 4 \pi^2$ at $k = \pi$. The mode has no horizontal velocity at $x = 0$ and $x = \pi / k$,
so a box of length $L = \pi H / k$ with impermeable side walls contains exactly one convection cell.

**Setup**

The domain is the box $[0, L] \times [0, H]$ with $L = H$, i.e. $k = \pi$.

| Boundary | Flow | Transport |
|----------|------|-----------|
| top, $z = H$ | no flow | $X = 1$ |
| bottom, $z = 0$ | no flow, pressure fixed at $(0, 0)$ | $X = 0$ |
| sides, $x = 0$ and $x = L$ | no flow | no flux |

The initial condition is $X = z/H + A_0 \cos(\pi x / L) \sin(\pi z / H)$ with $A_0 = 10^{-3}$ and the
hydrostatic pressure of the linear profile.

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Height, length | $H$, $L$ | 1, 1 | m |
| Reference density | $\varrho$ | 1000 | kg/m³ |
| Density increase at $X = 1$ | $\Delta\varrho$ | 500 | kg/m³ |
| Viscosity | $\mu$ | $10^{-3}$ | Pa s |
| Molecular diffusion coefficient | $D$ | $10^{-9}$ | m²/s |
| Porosity | $\phi$ | 0.4 | - |
| Tortuosity | $\tau$ | 1 | - |
| Permeability | $K$ | $6.44 \times 10^{-15}$ | m² |
| Gravitational acceleration | $g$ | 9.81 | m/s² |

The derived parameters are $D_\mathrm{eff} = 4 \times 10^{-10}$ m²/s, the time unit
$\phi H^2 / D_\mathrm{eff} = 10^9$ s and $Ra = 8 \pi^2 = 2 Ra_c$, for which $\sigma = 2 \pi^2$.
Since $\sigma$ depends on the parameters only through $Ra$, the density contrast $\Delta\varrho / \varrho = 0.5$
does not change the growth rate of the Boussinesq model. It makes the result sensitive to the density in
the storage and advection terms, which has to stay at the reference density $\varrho$.

**Discretization and measurement of the growth rate**

The @ref OnePNCModel with mass fractions is discretized with the box scheme (@ref BoxDiscretization) on a
structured grid. The total mass balance replaces the balance of the solvent, so that the boundaries with
fixed mass fraction are impermeable. The constant-tortuosity model gives $D_\mathrm{eff} = \phi \tau D$.
Time integration uses the implicit Euler method with 20 steps of size $\Delta t = 10^{-3}$ time units,
Newton's method and a direct linear solver.

After each time step, the amplitude $A_n$ is the projection of $X - z/H$ onto the mode
$\cos(\pi x / L) \sin(\pi z / H)$, evaluated with the sub-control volumes as quadrature. An implicit Euler
step multiplies a mode with growth rate $\sigma_h$ by $1 / (1 - \sigma_h \Delta t)$, so

$$\sigma_h = \frac{1 - A_{n-1} / A_n}{\Delta t}$$

is the growth rate of the spatially discretized problem, without time discretization error. The error
$|\sigma_h - \sigma| / (k^2 + \pi^2)$ is measured at the last time step, relative to the diffusive decay
rate $k^2 + \pi^2$, which does not vanish at $Ra = Ra_c$.

**Results**

The figures and the tables below are produced by running

```bash
python3 run_benchmark.py ./test_1p2c_hortonrogerslapwood_box
```

in the build directory `test/porousmediumflow/1pnc/1p2c/isothermal/hortonrogerslapwood`. The images of this
page are regenerated by passing `doc/doxygen/images` of the DuMux source directory as second argument.
The script varies the Rayleigh number through the permeability.

![Perturbation of the linear profile and streamlines of the Darcy velocity](hortonrogerslapwood_cells.png)

The figure above shows the perturbation $X - z/H$ of the linear profile, normalized by its maximum, and
streamlines of the Darcy velocity after 20 time steps at $Ra = 2 Ra_c$, for $k = \pi$ ($64 \times 64$ cells)
and for $k = \pi/2$ in the box $[0, 2H] \times [0, H]$ ($128 \times 64$ cells). Fluid with excess solute (red)
is denser and sinks, depleted fluid (blue) rises. The box of length $2H$ contains a single elongated cell:
the initial condition contains only the mode with $k = \pi/2$, and the mode with $k = \pi$, which forms two
cells and grows faster at this Rayleigh number ($\sigma = 2\pi^2$ against $\sigma = 3.45$), stays at round-off
level (relative amplitude $6 \times 10^{-13}$ after 20 time steps). The normalized perturbation differs from the mode shape
$\cos(k x) \sin(\pi z / H)$ by at most $1.1 \times 10^{-2}$ for $k = \pi$ and $4.8 \times 10^{-3}$ for $k = \pi/2$.
For $k = \pi$, the difference is $2.1 \times 10^{-2}$, $1.1 \times 10^{-2}$ and $5.9 \times 10^{-3}$ on
$32 \times 32$, $64 \times 64$ and $128 \times 128$ cells, and $2.0 \times 10^{-2}$, $1.0 \times 10^{-2}$ and
$5.0 \times 10^{-3}$ with the smaller initial amplitude $A_0 = 10^{-4}$, so it is a discretization error.
It comes from the upwind weighting of the advective flux: with central weighting (`Flux.UpwindWeight = 0.5`),
the difference for $k = \pi/2$ drops from $4.8 \times 10^{-3}$ to $3.4 \times 10^{-4}$, while the growth rate changes
from 3.4451 to 3.4463 (exact 3.4544).

![Growth rate over the Rayleigh number (left) and error under grid refinement (right)](hortonrogerslapwood_growthrate.png)

For $Ra = 8 \pi^2$ and $k = \pi$, the growth rate converges with second order under grid refinement:

| Grid | $\sigma_h$ | error | rate |
|------|------------|-------|------|
| $8 \times 8$ | 18.02445 | $8.69 \times 10^{-2}$ | - |
| $16 \times 16$ | 19.28944 | $2.28 \times 10^{-2}$ | 1.93 |
| $32 \times 32$ | 19.62142 | $5.97 \times 10^{-3}$ | 1.93 |
| $64 \times 64$ | 19.70738 | $1.61 \times 10^{-3}$ | 1.89 |
| $128 \times 128$ | 19.73002 | $4.66 \times 10^{-4}$ | 1.79 |

The exact growth rate is $\sigma = 2 \pi^2 = 19.73921$. For $0.5 \leq Ra / Ra_c \leq 4$ on $64 \times 64$
cells, the error is between $2.0 \times 10^{-5}$ and $4.7 \times 10^{-3}$. In the box $[0, 2H] \times [0, H]$
with $128 \times 64$ cells, the mode with $k = \pi/2$ grows for $Ra / Ra_c > 1.5625$, and the error is
between $4.2 \times 10^{-5}$ and $1.9 \times 10^{-3}$. Without perturbation ($A_0 = 0$), the maximum deviation
from the linear profile after 20 time steps is $2.2 \times 10^{-16}$: the discrete solution at rest is a
steady state.

The test `test_1p2c_hortonrogerslapwood_box` runs the setup above on $32 \times 32$ cells and fails if the
error exceeds $10^{-2}$. The test `test_1p2c_hortonrogerslapwood_decay_box` uses $Ra = 2 \pi^2$, for which
the perturbation decays with $\sigma = -\pi^2$; its error is $1.3 \times 10^{-4}$.
