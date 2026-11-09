# Onset of One-Sided Convection {#benchmark-one-sided-convection}

## Growth of a perturbation of a diffusive boundary layer with the Boussinesq approximation

**Problem description**

A porous layer of height $H$ between impermeable walls is saturated at its top boundary with a solute that makes
the fluid denser, for example CO<sub>2</sub> dissolving from a gas cap into water. The solute enters by diffusion
and forms a dense boundary layer whose thickness grows with the square root of time. At first, small perturbations
of the boundary layer decay. Once the layer is thick enough, they grow, which is the onset of convection
@cite Riaz2006 @cite Slim2014.

The benchmark covers this linear stage. A single perturbation mode with the wavenumber $k$ is added to the diffusive
boundary layer, and its profile and growth are compared with the linearized problem of Riaz et al. (2006)
@cite Riaz2006, solved independently of DuMux. The benchmark tests Darcy's law with the Boussinesq approximation,
`CVFEBoussinesqDarcysLaw`, where the density varies only in the gravity term,

$$\nabla \cdot \mathbf{v} = 0,
\qquad
\mathbf{v} = -\frac{K}{\mu} \left( \nabla p - (\varrho + \Delta\varrho X) \mathbf{g} \right),
\qquad
\phi \frac{\partial X}{\partial t} + \nabla \cdot \left( X \mathbf{v} - D_\mathrm{eff} \nabla X \right) = 0,$$

with the Darcy velocity $\mathbf{v}$, pressure $p$, solute mass fraction $X$, permeability $K$, viscosity $\mu$,
reference density $\varrho$, density increase $\Delta\varrho$ at $X = 1$, porosity $\phi$ and effective diffusion
coefficient $D_\mathrm{eff} = \phi \tau D$.

**Scaling and linearized problem**

Lengths are measured in units of $H$, velocities in units of $U = K \Delta\varrho g / \mu$ and time in units of
$\phi H / U$. The only parameter is the Rayleigh number $Ra = U H / D_\mathrm{eff}$. Let $t$ be the dimensionless time
since the top boundary is set to $X = 1$ and $\eta$ the dimensionless depth below the top. The base state is the
diffusive profile $C_0 = \mathrm{erfc}\left(\eta / (2 \sqrt{t / Ra})\right)$. A perturbation
$X - C_0 = \hat c(\eta, t) \cos(k x)$ causes the downward velocity $\hat w(\eta, t) \cos(k x)$ with @cite Riaz2006

$$\left( \frac{\partial^2}{\partial \eta^2} - k^2 \right) \hat w = -k^2 \hat c,
\qquad
\frac{\partial \hat c}{\partial t} = \frac{1}{Ra} \left( \frac{\partial^2}{\partial \eta^2} - k^2 \right) \hat c
- \frac{\partial C_0}{\partial \eta} \hat w,$$

and $\hat w = 0$ at $\eta = 0, 1$, $\hat c = 0$ at $\eta = 0$, $\partial \hat c / \partial \eta = 0$ at $\eta = 1$.
Since $C_0$ changes with time, there is no single growth rate (unlike in @ref benchmark-horton-rogers-lapwood): the
perturbation first decays, then grows, and depends on its initial shape. We use the shape of Riaz et al. (2006),
$\hat c(\eta, t_0) = \varepsilon\, \zeta \exp(-\zeta^2)$ with $\zeta = \eta / (2 \sqrt{t_0 / Ra})$. The results do not
depend on $\varepsilon$ (linear problem) or $t_0$. We compare the profile $\hat c(\eta, t)$, the growth rate
$\sigma = \mathrm{d} \ln \max_\eta |\hat c| / \mathrm{d}t$ plotted by Riaz et al., and the amplification
$A(t) / A(t_0)$ of the norm $A = \left( \int_0^1 \hat c^2 \,\mathrm{d}\eta \right)^{1/2}$, which is used by the test.

**Reference solution**

The script `linear_stability.py` solves the linearized problem independently of DuMux. The boundary layer stays thin
at these times, with a thickness $\sqrt{t / Ra} \leq 0.09\, H$ for $t \leq 4$, so the bottom has a negligible influence
and the layer can be treated as semi-infinite. There, $\hat c$ and $\hat w$ are expanded in the Laguerre functions
$\varphi_m(\eta) = e^{-\xi/2} \left( L_{m+1}(\xi) - L_m(\xi) \right)$, $\xi = \eta / s$, with a length scale $s$.
Since $L_n(0) = 1$, the functions vanish at the top, and they decay exponentially with depth, so they satisfy the
boundary conditions by construction. In time, the system is integrated with the Radau IIA method.
80, 100 and 120 functions give $\ln\left(A(2) / A(t_0)\right) = 1.35860074$ in all digits shown.

**Setup**

The domain is the box $[0, L] \times [0, H]$ with $L = H / 8$, i.e. $k = \pi H / L = 8 \pi$, one of the cases of
Figure 4 of Riaz et al. (2006) @cite Riaz2006. The box contains half a wavelength, and the impermeable side walls are
compatible with the mode. The top boundary has $X = 1$, all other boundaries are closed (pressure fixed at the bottom
left). The initial condition at $t_0 = 0.01$ is $X = C_0(\eta, t_0) + \varepsilon \cos(k x)\, \zeta \exp(-\zeta^2)$
with $\varepsilon = 10^{-6}$ and the hydrostatic pressure of $C_0$.

| Parameter | Value | Parameter | Value |
|-----------|-------|-----------|-------|
| $H$, $L$ | 1 m, 0.125 m | $\phi$, $\tau$ | 0.4, 1 |
| $\varrho$, $\Delta\varrho$ | 1000, 500 kg/m³ | $K$ | $4.077 \times 10^{-14}$ m² |
| $\mu$ | $10^{-3}$ Pa s | $D$ | $10^{-9}$ m²/s |

This gives $D_\mathrm{eff} = 4 \times 10^{-10}$ m²/s, $U = 2.0 \times 10^{-7}$ m/s and $Ra = 500$.

The @ref OnePNCModel is discretized with the box scheme (@ref BoxDiscretization) on a tensor-product grid with
$N$ cells in $x$ and $N$ cells in $0 \leq z \leq 0.7 H$ and $4N$ cells in $0.7 H \leq z \leq H$, graded towards the top.
The advective flux uses central differences (`Flux.UpwindWeight = 0.5`). Time integration uses the implicit Euler
method. After each time step, the amplitude $\hat c_h(\eta_j)$ on each row of vertices is the projection of $X$ minus its
row mean onto $\cos(k x)$. From it follow the norm $A_h$ and the maximum $M_h = \max_j |\hat c_h(\eta_j)|$, whose growth
rate $\sigma_h$ is computed by finite differences in time.

**Results**

The figures and tables are produced by `python3 run_benchmark.py ./test_1p2c_onesidedconvection_box` in the build
directory (pass `doc/doxygen/images` of the source directory as second argument to regenerate the images). The runs
use $16 \times 80$ cells and the time step $2 \times 10^{-3}$.

![Figure 1: profile of the perturbation over the depth at five times, linearized problem and DuMux](onesidedconvection_profile.png)

Figure 1 shows the normalized profile $\hat c(\eta, t) / \max_\eta |\hat c|$ at $t = 0.1$ to $2$. The maximum moves
down with the boundary layer. The largest difference of the normalized profiles between DuMux and the linearized
problem is $1.7 \times 10^{-2}$ at $t = 0.1$, when the boundary layer spans only a few cells, and below
$10^{-3}$ for $t \geq 1$.

![Figure 2: growth rate of the maximum of the perturbation for two wavenumbers, linearized problem, Riaz et al. (2006) and DuMux](onesidedconvection_growthrate.png)

Figure 2 shows $\sigma$ for $k = 4\pi$ and $8\pi$ at $Ra = 500$. The lines are the linearized problem, the triangles are
the markers of the initial-value solution in Figure 4 of Riaz et al. @cite Riaz2006 (`riaz2006_fig4_ivp.csv`), the dots are
DuMux in boxes of the width $L = H/4$ and $H/8$.

| $t$ | 1 | 2 | 4 |
|------|------|------|------|
| $k = 4\pi$, linearized problem / DuMux | 0.820 / 0.813 | 1.356 / 1.351 | 1.537 / 1.532 |
| $k = 8\pi$, linearized problem / DuMux | 1.849 / 1.838 | 1.981 / 1.972 | 1.740 / 1.732 |

For $t > 0.8$, DuMux differs from the linearized problem by at most 0.028, and the markers of Riaz et al. by at most
0.055. The latter is mainly a shift in time of 0.015 to 0.025 whose cause is not known.

The amplification $\ln (A(t) / A(t_0))$ of the linearized problem is $-1.630$, $-0.723$ and $1.359$ at $t = 0.5$, 1 and 2.
The time step is refined with the square of the cell size, so that the first-order error in time decreases like the
second-order error in space. The error of $\ln A$ at $t = 2$ is $-0.0846$, $-0.0213$ and $-0.0053$ on the grids
$8 \times 40$, $16 \times 80$ and $32 \times 160$ with $\Delta t = 4 \times 10^{-3}$, $10^{-3}$ and
$2.5 \times 10^{-4}$, which is a convergence rate of 2.0.

**Test**

The test `test_1p2c_onesidedconvection_box` runs the setup with $N = 16$ and $\Delta t = 2 \times 10^{-3}$ and fails if
$|\ln (A_h(2) / A_h(t_0)) - 1.35860074|$ exceeds 0.025; the error is 0.018.
