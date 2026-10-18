# Low-Rayleigh Elder Problem {#benchmark-elder}

## Solute plume below a source on top of a porous layer with the Boussinesq approximation

**Problem description**

Elder @cite Elder1967 heated part of the base of a Hele-Shaw cell and observed convection plumes.
Voss and Souza @cite VossSouza1987 formulated the problem for solute transport: a closed porous box of height $H$
and width $4H$, initially free of solute, takes up dense brine through the central half of its top boundary,
and the solute concentration at the bottom is zero. With the original parameters, the Rayleigh number defined
below is ${Ra = 400}$, and the problem is a common benchmark for density-driven flow @cite DierschKolditz2002.
At ${Ra = 400}$, the problem has at least three stable steady states, with one, two and three downward plumes
@cite Johannsen2003 @cite vanReeuwijk2009. Which of them a simulation reaches depends on the grid, the time step
and the scheme, so the original problem does not measure discretization errors. With a pseudo-spectral method,
van Reeuwijk et al. @cite vanReeuwijk2009 found that a single steady state exists for ${Ra < 76}$, and proposed
the problem at ${Ra = 60}$ as a benchmark, the low-Rayleigh Elder problem. Here, DuMux is compared with a Python
port of their pseudo-spectral solver @cite vanReeuwijkMathias2019. The benchmark tests `CVFEBoussinesqDarcysLaw`
with the box scheme in a transient problem in which the flow depends on the solute distribution.

With the Boussinesq approximation, the density varies only in the gravity term of Darcy's law,

$$\nabla \cdot \mathbf{v} = 0,
\qquad
\mathbf{v} = -\frac{K}{\mu} \left( \nabla p - (\varrho + \Delta\varrho X) \mathbf{g} \right),
\qquad
\phi \frac{\partial X}{\partial t} + \nabla \cdot \left( X \mathbf{v} - D_\mathrm{eff} \nabla X \right) = 0,$$

where $\mathbf{v}$ is the Darcy velocity, $p$ the pressure, $X$ the solute mass fraction, $K$ the permeability,
$\mu$ the viscosity, $\varrho$ the constant reference density of the fluid, $\Delta\varrho$ the density increase
at ${X = 1}$, $\mathbf{g}$ the gravitational acceleration with magnitude $g$, $\phi$ the porosity and
${D_\mathrm{eff} = \phi \tau D}$ the effective diffusion coefficient with tortuosity $\tau$ and molecular
diffusion coefficient $D$. The problem returns the density deviation $\Delta\varrho X$ with
`buoyantDensity(volVars, phaseIdx)`. With lengths in units of $H$ and time in units of
${\phi H^2 / D_\mathrm{eff}}$, the only parameter is the Rayleigh number

$$Ra = \frac{K \Delta\varrho g H}{\mu D_\mathrm{eff}}.$$

**Setup**

The problem is symmetric about the vertical center line of the box. As in the reference solution, only the
right half ${[0, L] \times [0, H]}$ with ${L = 2H}$ is computed, and the center line ${x = 0}$ is a boundary
without fluid and solute flux.

| Boundary | Flow | Transport |
|----------|------|-----------|
| top, ${z = H}$ | no flow, pressure fixed at ${(L, H)}$ | ${X = f(x)}$ |
| bottom, ${z = 0}$ | no flow | ${X = 0}$ |
| center line ${x = 0}$, side ${x = L}$ | no flow | no flux |

In the original problem, $f$ is a step from 1 to 0 at the edge ${x_s = H}$ of the source. The step makes the
solute flux at ${x = x_s}$ infinite, and no discretization converges pointwise near it. Following
van Reeuwijk et al., the step is replaced by

$$f(x) = \frac{1}{2} \left( 1 - \tanh\left( \frac{x - x_s}{\delta} \right) \right)$$

with the transition length $\delta$. They use ${\delta = 0.05 H}$ for their results at ${Ra = 400}$ and
${\delta = 0.005 H}$ at ${Ra = 60}$. This benchmark uses ${\delta = 0.05 H}$, for which the reference solution
converges pointwise; results for ${\delta = 0.005 H}$ and the step are given below. The initial condition is
${X = 0}$ with the hydrostatic pressure ${p = p_0 + \varrho g (H - z)}$, ${p_0 = 10^5}$ Pa.

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Height, length of the half domain | $H$, $L$ | 150, 300 | m |
| Half-width of the source | $x_s$ | 150 | m |
| Transition length of the source | $\delta$ | 7.5 | m |
| Reference density | $\varrho$ | 1000 | kg/m³ |
| Density increase at ${X = 1}$ | $\Delta\varrho$ | 30 | kg/m³ |
| Viscosity | $\mu$ | ${10^{-3}}$ | Pa s |
| Molecular diffusion coefficient | $D$ | ${3.565 \times 10^{-6}}$ | m²/s |
| Porosity | $\phi$ | 0.1 | - |
| Tortuosity | $\tau$ | 1 | - |
| Permeability | $K$ | ${4.8454 \times 10^{-13}}$ | m² |
| Gravitational acceleration | $g$ | 9.81 | m/s² |

The derived parameters are ${D_\mathrm{eff} = 3.565 \times 10^{-7}}$ m²/s, the time unit
${\phi H^2 / D_\mathrm{eff} = 6.311 \times 10^9}$ s, or 200 years, and ${Ra = 60}$. The porosity, viscosity, density,
diffusion coefficient and gravitational acceleration are those of Voss and Souza @cite VossSouza1987, whose
diffusive flux $\phi \varrho D \nabla X$ corresponds to ${\tau = 1}$. The permeability gives
${Ra = 2 \Delta\varrho / (\mathrm{kg/m^3})}$ exactly, so the density increase ${\Delta\varrho = 200}$ kg/m³ of the
original problem gives ${Ra = 400}$; their value ${K = 4.845 \times 10^{-13}}$ m² gives ${Ra = 399.97}$.

**Reference solution**

van Reeuwijk et al. @cite vanReeuwijk2009 write the dimensionless equations with the stream function $\psi$,

$$\frac{\partial c}{\partial t} + Ra \left( \frac{\partial (q_x c)}{\partial x} + \frac{\partial (q_z c)}{\partial z} \right)
= \nabla^2 c,
\qquad
\nabla^2 \psi = \frac{\partial c}{\partial x},
\qquad
q_x = \frac{\partial \psi}{\partial z}, \quad q_z = -\frac{\partial \psi}{\partial x},$$

where ${c = X}$, and the dimensionless Darcy velocity ${(q_x, q_z)}$ is ${\mathbf{v} H / (Ra D_\mathrm{eff})}$.
The concentration is expanded in a cosine series and the stream function in a sine series in $x$, and both in
Chebyshev polynomials in $z$. The transport equation is collocated at the Chebyshev-Gauss-Lobatto points and the
stream function of each horizontal mode is computed with the Chebyshev tau method. `elder_pseudospectral.py`
ports this method to Python. It solves the problem of the paper; the MATLAB code of the authors
@cite vanReeuwijkMathias2019 runs a demonstration on the domain ${[0, \pi] \times [-1, 1]}$ with a different
source and initial condition. Time integration uses the implicit-explicit Runge-Kutta scheme ARS(4,4,3)
@cite AscherRuuthSpiteri1997 with implicit diffusion and explicit advection.

The port is verified as follows; the distances are in units of $H$.

- Without buoyancy (${Ra = 0}$), each mode $f_m \cos(k_m x)$ of the cosine series of $f$, with ${k_m = m \pi H / L}$,
  diffuses into the layer independently, and

  $$c(x, z, t) = \sum_{m \geq 0} f_m \cos(k_m x) \left( \frac{\sinh(k_m z)}{\sinh(k_m)}
  - \sum_{n \geq 1} \frac{2 (-1)^{n+1} n \pi}{k_m^2 + n^2 \pi^2} \sin(n \pi z) e^{-(k_m^2 + n^2 \pi^2) t} \right),$$

  where the first term is $z$ for ${m = 0}$. For the step, this agrees with the solution of van Reeuwijk et al. in
  Laplace space, inverted numerically, to ${1.9 \times 10^{-11}}$. With ${385 \times 49}$ collocation points and a time
  step of ${2.5 \times 10^{-4}}$, the port reproduces it with a maximum error of ${1.1 \times 10^{-12}}$ at ${t = 1}$ and
  ${2.1 \times 10^{-8}}$ at ${t = 0.03}$, where the time step limits the accuracy.
- At ${Ra = 60}$ with ${\delta = 0.005 H}$ and ${97 \times 33}$ points, as in the paper, the isolines
  ${c = 0.1, 0.2, \dots, 0.9}$ of the port lie within a median distance of ${4.8 \times 10^{-4}}$ to
  ${8.2 \times 10^{-4}}$ and a maximum distance of ${3.0 \times 10^{-3}}$ of the published isolines at the four times
  of the paper, read from the vector graphics of the article. The published isolines without buoyancy deviate
  from the exact solution by a median ${6.1 \times 10^{-4}}$ to ${7.8 \times 10^{-4}}$ and at most
  ${1.1 \times 10^{-2}}$, which is the precision of this reading.
- At ${Ra = 400}$ with ${\delta = 0.05 H}$, the port reaches the steady states with one and two downward plumes
  from the initial conditions given in the paper, with ${Sh = 3.4670}$ and ${Sh = 3.9785}$. The paper reports 3.5
  and 4.0, and its bifurcation diagram 3.4695 and 3.9789. The state with two plumes ceases to exist between
  ${Ra = 74}$ and ${Ra = 78}$; the paper finds ${Ra = 76}$.
- With the resolution used here, ${129 \times 41}$ points and a time step of ${2.5 \times 10^{-5}}$, the solute mass
  differs from that with ${193 \times 65}$ points and half the time step by less than ${10^{-10}}$, the Sherwood number
  by less than ${5 \times 10^{-7}}$, and the relative $L^2$ difference of $c$ is at most ${7.2 \times 10^{-9}}$.

**Discretization and measured quantities**

The @ref OnePNCModel with mass fractions is discretized with the box scheme (@ref BoxDiscretization) on a
structured grid with ${2N \times N}$ cells. The total mass balance replaces the balance of the solvent, so that
the boundaries with fixed mass fraction are impermeable. Time integration uses the three-stage diagonally
implicit Runge-Kutta scheme of Alexander @cite Alexander1977 (`Experimental::MultiStage::DIRKThirdOrderAlexander`)
with a time step of ${1.25 \times 10^{-3}}$ time units, Newton's method and a direct linear solver. The advective
flux is weighted with full upwinding (`Flux.UpwindWeight = 1`, the default) or centrally
(`Flux.UpwindWeight = 0.5`).

At the times ${t = 0.03}$, 0.1, 0.3 and 1 of the reference, which are 6, 20, 60 and 200 years, the run reports

- the solute mass ${M = \frac{1}{H^2} \int X \, \mathrm{d}A}$,
- the Sherwood number ${Sh = \frac{H}{L} \int_0^L \frac{\partial X}{\partial z} \Big|_{z = H} \mathrm{d}x}$,
  the solute flux through the top relative to the diffusive flux ${\varrho D_\mathrm{eff} L / H}$, and the
  corresponding $Sh_b$ at the bottom,
- the relative error ${\varepsilon = \| X_h - X_\mathrm{ref} \| / \| X_\mathrm{ref} \|}$ of the DuMux solution $X_h$
  against the reference solution $X_\mathrm{ref}$, where ${\| \cdot \|}$ is the $L^2$ norm over the domain, evaluated
  with the vertex values and the sub-control volumes as quadrature. (van Reeuwijk et al. measure the error in the
  $L^1$ norm.)

The mass fraction at the top and the bottom is fixed, so the solute mass in the control volumes of these vertices
is constant, and the fluxes through the top and the bottom are the sums of the fluxes over the faces of these
control volumes.

**Results**

The figures and tables below are produced by running

```bash
python3 run_benchmark.py ./test_1p2c_elder_box
```

in the build directory `test/porousmediumflow/1pnc/1p2c/isothermal/elder`. The images of this page are
regenerated by passing `doc/doxygen/images` of the DuMux source directory as second argument. The runs take
20 CPU hours, of which the grids with ${512 \times 256}$ cells take about 3 hours each; `--jobs N` sets the number
of parallel runs (8 by default) and `--max-cells 256` leaves out the finest grid. The script also computes the
pseudo-spectral reference solutions, stores every run with the commit and a hash of the uncommitted changes, and
reuses runs of the same executable and arguments.

![Solute mass fraction of DuMux and isolines of DuMux and of the reference solution](elder_concentration.png)

The figure above shows the solute mass fraction on ${256 \times 128}$ cells with central weighting at the four times
and the isolines ${X = 0.1, 0.2, \dots, 0.9}$ of DuMux (white) and of the reference solution (black, dashed). The relative $L^2$ error of this solution is at most ${1.7 \times 10^{-4}}$.

With central weighting, the error converges with approximately second order towards the reference solution at all
four times (observed orders in parentheses):

| Grid | $\varepsilon$, ${t = 0.03}$ | $\varepsilon$, ${t = 0.1}$ | $\varepsilon$, ${t = 0.3}$ | $\varepsilon$, ${t = 1}$ |
|------|------|------|------|------|
| ${64 \times 32}$ | ${2.26 \times 10^{-3}}$ | ${2.64 \times 10^{-3}}$ | ${1.69 \times 10^{-3}}$ | ${1.96 \times 10^{-3}}$ |
| ${128 \times 64}$ | ${5.84 \times 10^{-4}}$ (1.95) | ${6.74 \times 10^{-4}}$ (1.97) | ${4.35 \times 10^{-4}}$ (1.96) | ${4.89 \times 10^{-4}}$ (2.00) |
| ${256 \times 128}$ | ${1.48 \times 10^{-4}}$ (1.98) | ${1.69 \times 10^{-4}}$ (2.00) | ${1.09 \times 10^{-4}}$ (2.00) | ${1.22 \times 10^{-4}}$ (2.00) |
| ${512 \times 256}$ | ${3.82 \times 10^{-5}}$ (1.95) | ${4.21 \times 10^{-5}}$ (2.01) | ${2.68 \times 10^{-5}}$ (2.03) | ${3.05 \times 10^{-5}}$ (2.00) |

With full upwinding, the default of DuMux, the error converges with approximately first order, and on
${512 \times 256}$ cells it is 26 to 145 times larger than with central weighting:

| Grid | $\varepsilon$, ${t = 0.03}$ | $\varepsilon$, ${t = 0.1}$ | $\varepsilon$, ${t = 0.3}$ | $\varepsilon$, ${t = 1}$ |
|------|------|------|------|------|
| ${64 \times 32}$ | ${9.06 \times 10^{-3}}$ | ${1.36 \times 10^{-2}}$ | ${2.81 \times 10^{-2}}$ | ${1.57 \times 10^{-2}}$ |
| ${128 \times 64}$ | ${4.23 \times 10^{-3}}$ (1.10) | ${6.64 \times 10^{-3}}$ (1.03) | ${1.50 \times 10^{-2}}$ (0.90) | ${7.72 \times 10^{-3}}$ (1.02) |
| ${256 \times 128}$ | ${2.04 \times 10^{-3}}$ (1.05) | ${3.28 \times 10^{-3}}$ (1.02) | ${7.70 \times 10^{-3}}$ (0.97) | ${3.84 \times 10^{-3}}$ (1.01) |
| ${512 \times 256}$ | ${1.00 \times 10^{-3}}$ (1.02) | ${1.63 \times 10^{-3}}$ (1.01) | ${3.88 \times 10^{-3}}$ (0.99) | ${1.91 \times 10^{-3}}$ (1.00) |

The reference values of the solute mass are ${M = 0.20565403}$, 0.41801512, 0.61388503 and 0.63434861, and those of
the Sherwood number ${Sh = 1.91136015}$, 1.52996376, 1.25285776 and 1.33086957 at ${t = 0.03}$, 0.1, 0.3 and 1.
The relative errors of DuMux with central weighting are

| Grid | ${t = 0.03}$ | ${t = 0.1}$ | ${t = 0.3}$ | ${t = 1}$ |
|------|------|------|------|------|
| ${64 \times 32}$ | ${3.3 \times 10^{-3}}$ | ${8.9 \times 10^{-4}}$ | ${9.9 \times 10^{-4}}$ | ${1.2 \times 10^{-3}}$ |
| ${128 \times 64}$ | ${8.2 \times 10^{-4}}$ | ${2.1 \times 10^{-4}}$ | ${2.4 \times 10^{-4}}$ | ${2.9 \times 10^{-4}}$ |
| ${256 \times 128}$ | ${2.1 \times 10^{-4}}$ | ${5.4 \times 10^{-5}}$ | ${6.0 \times 10^{-5}}$ | ${7.1 \times 10^{-5}}$ |
| ${512 \times 256}$ | ${5.7 \times 10^{-5}}$ | ${1.5 \times 10^{-5}}$ | ${1.5 \times 10^{-5}}$ | ${1.8 \times 10^{-5}}$ |

for the solute mass and

| Grid | ${t = 0.03}$ | ${t = 0.1}$ | ${t = 0.3}$ | ${t = 1}$ |
|------|------|------|------|------|
| ${64 \times 32}$ | ${4.4 \times 10^{-4}}$ | ${1.6 \times 10^{-3}}$ | ${4.6 \times 10^{-4}}$ | ${1.5 \times 10^{-3}}$ |
| ${128 \times 64}$ | ${1.3 \times 10^{-4}}$ | ${4.2 \times 10^{-4}}$ | ${1.3 \times 10^{-4}}$ | ${3.9 \times 10^{-4}}$ |
| ${256 \times 128}$ | ${4.1 \times 10^{-5}}$ | ${1.1 \times 10^{-4}}$ | ${3.4 \times 10^{-5}}$ | ${9.7 \times 10^{-5}}$ |
| ${512 \times 256}$ | ${1.8 \times 10^{-5}}$ | ${2.7 \times 10^{-5}}$ | ${9.0 \times 10^{-6}}$ | ${2.4 \times 10^{-5}}$ |

for the Sherwood number. The Sherwood number at ${t = 0.03}$ converges more slowly on the finest grid because of the
time step: on ${64 \times 32}$ cells, the relative time error at ${t = 0.03}$ is ${4.3 \times 10^{-6}}$ for $M$ and
${9.4 \times 10^{-6}}$ for $Sh$, and at most ${1.2 \times 10^{-6}}$ at the later times. With full upwinding on
${512 \times 256}$ cells, the relative errors are ${6.7 \times 10^{-4}}$ ($M$) and ${1.0 \times 10^{-3}}$ ($Sh$)
at ${t = 0.03}$ and ${8.3 \times 10^{-4}}$ and ${1.4 \times 10^{-3}}$ at ${t = 1}$.

![Error under grid refinement (left) and time error (right)](elder_convergence.png)

The time error is measured on ${64 \times 32}$ cells with central weighting as the relative $L^2$ difference to DIRK3
with a time step of ${3.125 \times 10^{-4}}$, maximized over the four times:

| $\Delta t$ | implicit Euler | DIRK3 |
|------|------|------|
| ${10^{-2}}$ | ${6.15 \times 10^{-2}}$ | ${1.66 \times 10^{-3}}$ |
| ${5 \times 10^{-3}}$ | ${3.22 \times 10^{-2}}$ (0.94) | ${2.51 \times 10^{-4}}$ (2.72) |
| ${2.5 \times 10^{-3}}$ | ${1.65 \times 10^{-2}}$ (0.97) | ${3.99 \times 10^{-5}}$ (2.65) |
| ${1.25 \times 10^{-3}}$ | ${8.34 \times 10^{-3}}$ (0.98) | ${6.16 \times 10^{-6}}$ (2.69) |
| ${6.25 \times 10^{-4}}$ | ${4.20 \times 10^{-3}}$ (0.99) |  |
| ${3.125 \times 10^{-4}}$ | ${2.11 \times 10^{-3}}$ (1.00) |  |

The implicit Euler method converges with approximately first order; at the time step of the grid study, its error
${8.3 \times 10^{-3}}$ is larger than the spatial error on ${64 \times 32}$ cells. DIRK3 converges with an order between
two and three, below its classical order 3, as expected for a method with stage order 1 applied to a stiff problem.
With implicit Euler, the solute mass changes in each step exactly by the time step times the fluxes through the top
and the bottom at the end of the step: the largest deviation is ${1.2 \times 10^{-14}}$ for ${\Delta t = 10^{-2}}$.

Without buoyancy (${\Delta\varrho = 0}$), the error against the exact solution given above converges with approximately
second order:

| Grid | $\varepsilon$, ${t = 0.03}$ | $\varepsilon$, ${t = 0.1}$ | $\varepsilon$, ${t = 0.3}$ | $\varepsilon$, ${t = 1}$ |
|------|------|------|------|------|
| ${64 \times 32}$ | ${8.95 \times 10^{-4}}$ | ${4.71 \times 10^{-4}}$ | ${3.87 \times 10^{-4}}$ | ${3.71 \times 10^{-4}}$ |
| ${128 \times 64}$ | ${2.23 \times 10^{-4}}$ (2.00) | ${1.18 \times 10^{-4}}$ (1.99) | ${9.72 \times 10^{-5}}$ (1.99) | ${9.32 \times 10^{-5}}$ (1.99) |
| ${256 \times 128}$ | ${5.45 \times 10^{-5}}$ (2.03) | ${2.95 \times 10^{-5}}$ (2.00) | ${2.42 \times 10^{-5}}$ (2.00) | ${2.33 \times 10^{-5}}$ (2.00) |
| ${512 \times 256}$ | ${1.27 \times 10^{-5}}$ (2.10) | ${7.34 \times 10^{-6}}$ (2.01) | ${6.05 \times 10^{-6}}$ (2.00) | ${5.81 \times 10^{-6}}$ (2.00) |

**Sharpness of the source**

For ${\delta = 0.005 H}$ and for the step, a reference solution with an error below that of the finest grids is not
available. For these, and for comparison ${\delta = 0.05 H}$, the differences of the mass fraction at ${t = 0.3}$ to the
next finer grid at the common vertices, with central weighting, are

| Grids | ${\delta = 0.05 H}$ | ${\delta = 0.005 H}$ | step |
|------|------|------|------|
| ${64 \times 32}$, ${128 \times 64}$ | ${1.26 \times 10^{-3}}$ | ${2.45 \times 10^{-3}}$ | ${2.47 \times 10^{-3}}$ |
| ${128 \times 64}$, ${256 \times 128}$ | ${3.26 \times 10^{-4}}$ (1.95) | ${7.79 \times 10^{-4}}$ (1.65) | ${9.41 \times 10^{-4}}$ (1.39) |
| ${256 \times 128}$, ${512 \times 256}$ | ${8.24 \times 10^{-5}}$ (1.98) | ${2.25 \times 10^{-4}}$ (1.79) | ${4.06 \times 10^{-4}}$ (1.21) |

in the relative $L^2$ norm and

| Grids | ${\delta = 0.05 H}$ | ${\delta = 0.005 H}$ | step |
|------|------|------|------|
| ${64 \times 32}$, ${128 \times 64}$ | ${4.80 \times 10^{-3}}$ | ${2.02 \times 10^{-2}}$ | ${2.03 \times 10^{-2}}$ |
| ${128 \times 64}$, ${256 \times 128}$ | ${1.37 \times 10^{-3}}$ (1.81) | ${1.68 \times 10^{-2}}$ (0.26) | ${1.97 \times 10^{-2}}$ (0.04) |
| ${256 \times 128}$, ${512 \times 256}$ | ${3.60 \times 10^{-4}}$ (1.93) | ${8.80 \times 10^{-3}}$ (0.93) | ${1.91 \times 10^{-2}}$ (0.04) |

in the maximum norm. For ${\delta = 0.05 H}$, both decrease with approximately second order. For the step, the largest
difference is at a vertex next to the edge of the source and stays at ${2 \times 10^{-2}}$ under refinement, so the
solution does not converge pointwise there. A difference of this size on an area of about one cell contributes a part
proportional to the cell size to the $L^2$ norm, and the order of the $L^2$ difference falls towards one. For
${\delta = 0.005 H}$, the maximum difference decreases from ${2.0 \times 10^{-2}}$ to ${1.7 \times 10^{-2}}$ and
${8.8 \times 10^{-3}}$; the finest grid is the first with cells smaller than $\delta$ (cell size ${H/256}$).

**Several steady states at Ra = 400**

With ${\Delta\varrho = 200}$ kg/m³, ${Ra = 400}$. The pseudo-spectral solution (${129 \times 41}$ points,
${\delta = 0.05 H}$) reaches a steady state with one downward plume on the center line (${Sh = 3.4670}$) from the
initial condition ${X = \frac{1}{2} (1 + \cos(\pi x / H))}$ for ${x < H}$ and ${X = 0}$ elsewhere, one of those given by
van Reeuwijk et al., and a steady state with a downward plume on each side of the center line (${Sh = 3.9785}$)
from ${X = 0}$.
DuMux with full upwinding, started from ${X = 0}$, reaches either of them at ${t = 1}$ depending on the grid and the
time step of DIRK3:

| Grid | ${\Delta t = {2.5 \times 10^{-3}}}$ | ${\Delta t = {1.25 \times 10^{-3}}}$ |
|------|------|------|
| ${64 \times 32}$ | two plumes, ${Sh = 4.0011}$ | one plume, ${Sh = 3.4581}$ |
| ${128 \times 64}$ | two plumes, ${Sh = 3.9908}$ | two plumes, ${Sh = 3.9908}$ |
| ${256 \times 128}$ | one plume, ${Sh = 3.4655}$ | two plumes, ${Sh = 3.9864}$ |

![Steady states at Ra = 400 of the reference solution (top) and of DuMux](elder_ra400.png)

The figure above shows the two steady states of the reference solution and one state of each kind reached by DuMux.
The Sherwood numbers of the DuMux states approach those of the reference states under refinement. Which state a
run reaches changes with the grid and with the time step, without a monotone pattern, so the runs at ${Ra = 400}$
do not form a sequence that converges to one solution.

**Tests**

The test `test_1p2c_elder_box` runs the setup above on ${64 \times 32}$ cells with central weighting and DIRK3 with a
time step of ${2.5 \times 10^{-3}}$ to ${t = 0.3}$ and fails if the relative error of $M$ or $Sh$ at ${t = 0.03}$, 0.1 or
0.3 against the reference values exceeds ${5 \times 10^{-3}}$; the largest error is ${3.3 \times 10^{-3}}$, of $M$ at
${t = 0.03}$. The test `test_1p2c_elder_nobuoyancy_box` runs the same setup without buoyancy and fails if the
relative $L^2$ error against the exact solution exceeds ${1.5 \times 10^{-3}}$; the largest error is
${8.8 \times 10^{-4}}$, at ${t = 0.03}$. Each test takes about 20 seconds.
