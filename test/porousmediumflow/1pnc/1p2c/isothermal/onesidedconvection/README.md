# Onset of One-Sided Convection in a Porous Layer {#benchmark-one-sided-convection}

## Diffusive boundary layer and growth of a finger mode with the Boussinesq approximation

**Problem description**

A porous layer of height $H$ between impermeable walls is saturated at its top boundary with a solute that makes
the fluid denser, for example CO<sub>2</sub> dissolving from a gas cap into water. The solute enters by diffusion
and forms a dense boundary layer whose thickness grows with the square root of time. Once the layer is thick
enough, fingers of dense fluid start to sink into the layer and enhance the dissolution
@cite Riaz2006 @cite Slim2014. The benchmark covers the stages before the fingers interact, which are
deterministic and have references independent of DuMux: the diffusive boundary layer without perturbation, compared
with the exact solution, and the growth of a single finger mode, compared with the linearized problem of Riaz et al.
(2006) @cite Riaz2006. It tests Darcy's law with the Boussinesq approximation, `CVFEBoussinesqDarcysLaw`, with a
boundary layer much thinner than the layer.

With the Boussinesq approximation, the density varies only in the gravity term of Darcy's law,

$$\nabla \cdot \mathbf{v} = 0,
\qquad
\mathbf{v} = -\frac{K}{\mu} \left( \nabla p - (\varrho + \Delta\varrho X) \mathbf{g} \right),
\qquad
\phi \frac{\partial X}{\partial t} + \nabla \cdot \left( X \mathbf{v} - D_\mathrm{eff} \nabla X \right) = 0,$$

where $\mathbf{v}$ is the Darcy velocity, $p$ the pressure, $X$ the solute mass fraction, $K$ the permeability,
$\mu$ the viscosity, $\varrho$ the constant reference density of the fluid, $\Delta\varrho$ the density increase at
$X = 1$, $\mathbf{g}$ the gravitational acceleration with magnitude $g$, $\phi$ the porosity and
$D_\mathrm{eff} = \phi \tau D$ the effective diffusion coefficient with tortuosity $\tau$ and molecular diffusion
coefficient $D$. The problem returns the density deviation $\Delta\varrho X$ with
`buoyantDensity(volVars, phaseIdx)`.

**Scaling**

Following Riaz et al. (2006) @cite Riaz2006, lengths are measured in units of $H$, velocities in units of the
buoyancy velocity $U = K \Delta\varrho g / \mu$ and time in units of $\phi H / U$. The only parameter of the problem
is the Rayleigh number

$$Ra = \frac{U H}{D_\mathrm{eff}} = \frac{K \Delta\varrho g H}{\mu D_\mathrm{eff}}.$$

In the following, $t$ is the dimensionless time and $\eta = d / H$ the dimensionless depth below the top boundary,
where $d = H - z$ for the vertical coordinate $z$.

**Diffusive boundary layer**

Without perturbation the fluid stays at rest, and $X$ is the diffusive profile $C_0(\eta, t)$. As long as the
boundary layer is thin compared with $H$, $C_0 = \mathrm{erfc}\left(\eta / (2 \sqrt{t / Ra})\right)$, and the
dissolution flux, the solute mass flux through the top boundary in units of $\varrho U$, is

$$F(t) = \frac{1}{Ra} \frac{\partial C_0}{\partial \eta}\Big|_{\eta = 0} = \frac{1}{\sqrt{\pi Ra\, t}}.$$

**Linearized problem**

A perturbation $X - C_0 = \hat c(\eta, t) \cos(k x)$ with wavenumber $k$ (in units of $1/H$) causes the downward
velocity $\hat w(\eta, t) \cos(k x)$ (in units of $U$). The linearized equations are @cite Riaz2006

$$\left( \frac{\partial^2}{\partial \eta^2} - k^2 \right) \hat w = -k^2 \hat c,
\qquad
\frac{\partial \hat c}{\partial t} = \frac{1}{Ra} \left( \frac{\partial^2}{\partial \eta^2} - k^2 \right) \hat c
- \frac{\partial C_0}{\partial \eta} \hat w,$$

with $\hat w = 0$ at $\eta = 0$ and $\eta = 1$, $\hat c = 0$ at $\eta = 0$ and $\partial \hat c / \partial \eta = 0$
at $\eta = 1$. Unlike the steady base state of the Horton-Rogers-Lapwood problem
(@ref benchmark-horton-rogers-lapwood), the base state $C_0$ changes with time. A perturbation therefore has no
single growth rate, and its growth depends on its initial shape. In the setup below, the perturbation first decays
and grows once the boundary layer has become unstable. The reference is the solution of this one-dimensional problem for the same initial
perturbation as in the DuMux setup, the dominant mode of Riaz et al. (2006),

$$\hat c(\eta, t_0) = \varepsilon\, \zeta \exp(-\zeta^2), \qquad \zeta = \frac{\eta}{2 \sqrt{t_0 / Ra}},$$

with amplitude $\varepsilon$. The benchmark compares the amplification $A(t) / A(t_0)$ of the norm
$A(t) = \left( \int_0^1 \hat c^2 \,\mathrm{d}\eta \right)^{1/2}$.

**Reference solution**

The script `linear_stability.py` solves the linearized problem with Chebyshev collocation in $\eta$ and the
Radau IIA method in time, with $C_0$ of the layer of finite depth as a series of image solutions. Running
`python3 linear_stability.py` prints the following checks:
- For $k = 0$, the perturbation $\eta\, t^{-3/2} \exp(-Ra\, \eta^2 / (4t))$ is an exact solution of the semi-infinite
  layer. With 128 Chebyshev points the computed perturbation at $t = 0.5$ deviates from it by
  $1.6 \times 10^{-9}$ relative to its maximum.
- With 128, 192 and 256 Chebyshev points, $\ln\left(A(2) / A(t_0)\right) = 1.35860074$ for the setup below, identical
  in all digits shown.
- For $Ra = 500$ the growth rate of $\max_\eta |\hat c|$ is 0.82, 1.36 and 1.54 at $t = 1, 2, 4$ for $k = 4\pi$, and
  1.85, 1.98 and 1.74 for $k = 8\pi$. These values agree with the initial-value solution in Figure 4 of Riaz et al.
  (2006) @cite Riaz2006 within the precision of reading the figure, about 0.02.

**Setup**

The domain is the box $[0, L] \times [0, H]$ with $L = H / 8$, i.e. $k = \pi H / L = 8 \pi$. Since
$\partial_x \cos(k x)$ vanishes at $x = 0$ and $x = L$, the impermeable side walls are compatible with the mode, and
the box contains half a wavelength: dense fluid sinks along the wall at $x = 0$ and rises along the wall at $x = L$.

| Boundary | Flow | Transport |
|----------|------|-----------|
| top, $z = H$ | no flow | $X = 1$ |
| bottom, $z = 0$ | no flow, pressure fixed at $(0, 0)$ | no flux |
| sides, $x = 0$ and $x = L$ | no flow | no flux |

The initial condition at $t_0 = 0.01$ is $X = C_0(\eta, t_0) + \varepsilon \cos(k x)\, \zeta \exp(-\zeta^2)$
with $\varepsilon = 10^{-6}$ and the hydrostatic pressure of $C_0$. At $t_0$ the boundary layer thickness
$2 \sqrt{t_0 / Ra}$ is $8.9 \times 10^{-3} H$.

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Height, length | $H$, $L$ | 1, 0.125 | m |
| Reference density | $\varrho$ | 1000 | kg/m³ |
| Density increase at $X = 1$ | $\Delta\varrho$ | 500 | kg/m³ |
| Viscosity | $\mu$ | $10^{-3}$ | Pa s |
| Molecular diffusion coefficient | $D$ | $10^{-9}$ | m²/s |
| Porosity | $\phi$ | 0.4 | - |
| Tortuosity | $\tau$ | 1 | - |
| Permeability | $K$ | $4.077 \times 10^{-14}$ | m² |
| Gravitational acceleration | $g$ | 9.81 | m/s² |

The derived parameters are $D_\mathrm{eff} = 4 \times 10^{-10}$ m²/s, $U = 2.0 \times 10^{-7}$ m/s, the time unit
$\phi H / U = 2.0 \times 10^{6}$ s and $Ra = 500$. The case $Ra = 500$, $k = 8\pi$ is one of the cases of
Figure 4 of Riaz et al. (2006) @cite Riaz2006.

**Discretization and measurement**

The @ref OnePNCModel with mass fractions is discretized with the box scheme (@ref BoxDiscretization). The total mass
balance replaces the balance of the solvent, so that the top boundary with fixed mass fraction is impermeable. The
grid is a tensor-product grid with $N$ cells in $x$, $N$ equal cells in $0 \leq z \leq 0.7 H$ and $4N$ cells in
$0.7 H \leq z \leq H$ that shrink towards the top by the factor $1.03^{16/N}$ per cell (smallest cell
$1.6 \times 10^{-3} H$ for $N = 16$). Time integration uses the implicit Euler method with a constant time step,
Newton's method and a direct linear solver.

After each time step, the amplitude $\hat c_h(\eta_j)$ on each row $j$ of vertices is the projection of $X$ onto
$\cos(k x)$, evaluated with the sub-control volumes as quadrature, and $A_h$ is its norm over the depth with the same
quadrature. The dissolution flux $F_h$ is the change of the mean solute mass fraction in the domain per time step.
On the box scheme this equals the solute flux out of the vertices with fixed mass fraction at the top. It is
compared with the mean of the exact flux over the time step.

**Results**

The tables and the figures below are produced by running

```bash
python3 run_benchmark.py ./test_1p2c_onesidedconvection_box
```

in the build directory `test/porousmediumflow/1pnc/1p2c/isothermal/onesidedconvection`. The images of this page are
regenerated by passing `doc/doxygen/images` of the DuMux source directory as second argument.

![Amplification of the finger mode and error of its logarithm under grid refinement](onesidedconvection_growth.png)

The mode first decays: the amplification of the reference has its minimum 0.159 at $t = 0.262$. Once the boundary
layer has become unstable, the mode grows, to $A(2) / A(t_0) = 3.89$. The reference values of
$\ln (A / A(t_0))$ are $-1.63016$, $-0.72291$ and $1.35860$ at $t = 0.5$, 1 and 2. With the time step
$\Delta t = 10^{-3}$, the error of $\ln A$ under grid refinement is

| Cells | $t = 0.5$ | $t = 1$ | $t = 2$ | $t = 2$, without time error | rate |
|------|------|------|------|------|------|
| $8 \times 40$ | -0.0117 | -0.0390 | -0.0921 | -0.0948 | - |
| $16 \times 80$ | -0.0027 | -0.0092 | -0.0211 | -0.0239 | 1.99 |
| $32 \times 160$ | -0.0004 | -0.0017 | -0.0032 | -0.0059 | 2.00 |

The column "without time error" subtracts the error of the implicit Euler method estimated from the time-step
refinement below, and the rate refers to this column. Before $t = 0.1$ the error has a peak of 0.031, 0.016 and
0.013 on $8 \times 40$, $16 \times 80$ and $32 \times 160$ cells. It comes from the time step, a tenth of $t_0$, and
from the boundary layer at $t_0$, which spans only a few cells.

On $16 \times 80$ cells, the error of $\ln A$ changes with the time step as

| Time step | $t = 0.5$ | $t = 1$ | $t = 2$ |
|------|------|------|------|
| $2 \times 10^{-3}$ | -0.0024 | -0.0083 | -0.0184 |
| $10^{-3}$ | -0.0027 | -0.0092 | -0.0211 |
| $5 \times 10^{-4}$ | -0.0029 | -0.0096 | -0.0225 |

The differences between successive time steps halve with the time step, as for a first-order method. The time
error at $\Delta t = 10^{-3}$ is 0.0003, 0.0008 and 0.0027 at $t = 0.5$, 1 and 2.

![Perturbation of the horizontal mean and streamlines of the Darcy velocity at t = 2 on 32 x 160 cells, and the amplitude profile compared with the linearized problem](onesidedconvection_mode.png)

At $t = 2$, dense fluid sinks along the wall at $x = 0$, with the largest amplitude $\hat c$ at the depth $0.06 H$.
The normalized amplitude profile $\hat c_h / \max |\hat c_h|$ differs from the reference by at most
$2.3 \times 10^{-3}$, $6.3 \times 10^{-4}$ and $2.0 \times 10^{-4}$ on $8 \times 40$, $16 \times 80$ and
$32 \times 160$ cells.

Without perturbation ($16 \times 80$ cells, $\Delta t = 10^{-3}$), the mass fraction in a row of vertices varies by
at most $1.9 \times 10^{-12}$ up to $t = 2$: the horizontally uniform state at rest is a solution of the discrete
problem, up to round-off. The relative error of the dissolution flux is $3.0 \times 10^{-2}$ after the first time
step, when the boundary layer spans only a few cells, at most $5.5 \times 10^{-4}$ for $t > 0.1$ and
$1.0 \times 10^{-4}$ at $t = 2$.

**Beyond the linear regime**

Later, the fingers interact, merge and form plumes, and the dynamics become chaotic: finger positions and the time at
which the first finger reaches the bottom differ between realizations of the initial perturbation. Slim (2014)
@cite Slim2014 reports a regime with a statistically constant dissolution flux $F = 0.017$ in the units used here, independent
of $Ra$.
Comparisons in that regime require ensembles of realizations or wide domains and are not part of this benchmark.

**Tests**

The test `test_1p2c_onesidedconvection_box` runs the setup above with $N = 16$ and $\Delta t = 2 \times 10^{-3}$ and
fails if $|\ln (A_h(2) / A_h(t_0)) - \ln (A(2) / A(t_0))|$ exceeds 0.025, with the reference value 1.35860074 from
`linear_stability.py`; the error is 0.018. The test `test_1p2c_onesidedconvection_diffusion_box` runs the setup
without perturbation until $t = 1$ and fails if the relative error of the dissolution flux at $t = 1$ exceeds
$10^{-3}$ or if the mass fraction in a row of vertices varies by more than $10^{-10}$; the values are
$2.7 \times 10^{-4}$ and $2.9 \times 10^{-13}$.
