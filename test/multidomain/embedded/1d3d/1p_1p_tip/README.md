# Vessel with a Closed Tip (1D-3D Embedded Coupling) {#benchmark-embedded-1d3d-vessel-tip}

## Convergence of the embedded coupling methods at a vessel end with an impermeable cap

**Problem description**

The model is that of @ref benchmark-embedded-1d3d-vessel, but the vessel ends inside the tissue. The
vessel center-line is $\Lambda = \{ (0, 0, s) : s \in [0, L] \}$ in the tissue block
$\Omega = [-1, 1] \times [-1, 1] \times [0, 1]$. The vessel enters the tissue through the face $z = 0$ and ends
at $z = L$. The end of the vessel is closed by an impermeable cap: no fluid leaves the vessel through its end,
and the exchange with the tissue takes place only through the lateral vessel surface. In the mixed-dimension
model this means

- zero axial flux at the tip, $\kappa_v \, \mathrm{d}p_v / \mathrm{d}s = 0$ at $s = L$,
- no source beyond the tip: the line source lies on $\Lambda$, the surface source on the lateral surface
  $r = R$, $0 \leq z \leq L$ (the cap disc at $z = L$ carries no source), and the kernel source in the cylinder
  $r \leq \varrho$, $0 \leq z \leq L$,

with $r = \sqrt{x^2 + y^2}$ the distance from the $z$-axis, $R$ the vessel radius and $\varrho$ the kernel
radius. The face $z = 0$ is a plane of symmetry (no flow), i.e. the setup is half of a vessel on
$s \in [-L, L]$ in an infinite tissue. Near the tip, the tissue pressure is no longer the two-dimensional
potential of the through-going vessel, and the three coupling methods (line source LS, cylinder surface source
CSS, distributed source DS) have different exact solutions outside the vessel as well. For the parameters below,
the exact tissue pressures of CSS and DS differ from that of LS by up to about $9.5 \times 10^{-3}$ and
$4.4 \times 10^{-3}$ at $r \geq R$ (3.3 % and 1.5 % of the local pressure), close to the vessel surface below
the tip.
Therefore each method is compared with its own exact solution.

![Setup in the plane y = 0 through the vessel axis with the boundary conditions (hatched: no flow). The exact solution is that of the vessel and its mirror image below the symmetry plane.](embedded1d3d_tip_setup.svg){html: width=85%}

**Exact solution**

The exchange rate per unit length is the constant $q = q_0$ and the vessel conductivity $\kappa_v$ is constant.
With zero axial flux at the tip, the vessel pressure is

$$p_v(s) = p_\mathrm{tip} + \frac{q_0}{2 \kappa_v} (L - s)^2 .$$

The tissue pressures are superpositions of the potential of a uniform source on a segment parallel to the
$z$-axis at distance $d$,

$$F(d, z) = \frac{q_0}{4\pi} \left( \operatorname{asinh} \frac{L - z}{d} + \operatorname{asinh} \frac{L + z}{d} \right),$$

which includes the mirror image of the vessel. The line source solution is $p_t^\mathrm{LS}(r, z) = F(r, z)$. The
surface and kernel solutions average $F$ over the circle of radius $R$ and over the disc of radius $\varrho$,

$$p_t^\mathrm{CSS}(r, z) = \frac{1}{\pi} \int_0^\pi F(d_R, z) \, \mathrm{d}\varphi,
\qquad
p_t^\mathrm{DS}(r, z) = \frac{1}{\pi} \int_0^\pi \int_0^\varrho \frac{2 r'}{\varrho^2} F(d_{r'}, z) \, \mathrm{d}r' \, \mathrm{d}\varphi,
\qquad
d_{a}^2 = r^2 + a^2 - 2 r a \cos\varphi .$$

The integral over $r'$ has a closed form. The remaining integral over $\varphi$, whose integrand is singular at
$\varphi = 0$ for points on the vessel surface or inside the kernel, is evaluated with tanh-sinh quadrature to a
relative tolerance of $10^{-13}$ (`exactsolution.hh`).

The exchange rate of each method is $q = \xi \beta(s) (p_v - \bar{p}_t)$, with the tissue pressure $\bar{p}_t$
of the method: the circle average at distance $R$ for LS and CSS, $\bar{p}_t = p_t(R, s)$, and the center-line
pressure for DS, $\bar{p}_t = p_t(0, s)$, which has a closed form. The flux scaling factor of DS is
$\xi = \left( 1 + \frac{\beta}{2\pi} \left( \ln \frac{\varrho}{R} - \frac12 \right) \right)^{-1}$
@cite Koch2020 (Eq. 2.12), and $\xi = 1$ for LS and CSS. The coefficient $\beta(s)$ is chosen such that

$$\xi \beta(s) \left( p_v(s) - \bar{p}_t(s) \right) = q_0,$$

so that $p_v$ and $q$ are the same for all methods. For DS, where $\xi$ depends on $\beta$, this is solved for
$\beta$ in closed form. For the parameters below, $\beta$ ranges from 1.276 to 1.503 (LS), 1.278 to 1.510 (CSS)
and 1.271 to 1.499 (DS) along the vessel. The discrete kernel method evaluates $\xi$ at the mean distance between
the vessel segment and the tissue cell @cite Koch2020 (Eqs. 3.18 and 3.20), which tends to the value above under
grid refinement. The kernel coupling manager passes the identifier of each coupling source to
`fluxScalingFactor(distance, radius, kernelWidth, id)` of the tissue problem, which evaluates $\beta$ at the
position of the source on the center-line.

| Boundary | Condition |
|----------|-----------|
| vessel, $s = 0$ | $p_v = p_\mathrm{tip} + q_0 L^2 / (2 \kappa_v)$ |
| vessel, $s = L$ | no flow |
| tissue, $z = 0$ | no flow |
| tissue, $z = 1$, $x = \pm 1$ and $y = \pm 1$ | exact tissue pressure |

**Parameters**

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Vessel length | $L$ | 0.6 |
| Vessel radius | $R$ | 0.1 |
| Kernel radius | $\varrho$ | $R$ (`MixedDimension.KernelWidthFactor = 1`) |
| Exchange rate per unit length | $q_0$ | 1 (`Problem.ExchangeRate`) |
| Vessel conductivity | $\kappa_v$ | 1 (`SpatialParams.VesselConductivity`) |
| Vessel pressure at the tip | $p_\mathrm{tip}$ | 1 (`Problem.TipPressure`) |
| Tissue permeability, fluid viscosity, fluid density | | 1 |

**Discretization**

The tissue grid has $c \times c \times c/2$ cells (`Dune::YaspGrid`) and the vessel grid $0.3 c$ cells
(`Dune::FoamGrid`), so that both have the cell size $h = 2/c$. The convergence study uses $c = 10, 20, 40, 80$,
i.e. $h/R = 2, 1, 0.5, 0.25$. On these grids the tip lies on a cell face and $R/h$ is an integer for $h \leq R$,
as in @ref benchmark-embedded-1d3d-vessel, where this alignment is needed for monotone convergence. The vessel
uses the cell-centered finite volume method with two-point flux approximation (@ref CCTpfaDiscretization), the
tissue the same method or the box method (@ref BoxDiscretization). The coupling is implemented by the coupling
managers of @ref EmbeddedCoupling with the coupling modes `Average` (LS), `Surface` (CSS) and `Kernel` (DS).

The errors are relative discrete $L^2$ errors with the exact solution evaluated at the degrees of freedom and the
control volumes $K$ as weights,

$$e_t = \left( \frac{\sum_K \lvert K \rvert \left( p_{t,K} - p_t(\mathbf{x}_K) \right)^2}{\sum_K \lvert K \rvert \, p_t(\mathbf{x}_K)^2} \right)^{1/2},
\qquad
e_q = \left( \frac{\sum_K \lvert K \rvert \left( \bar{q}_K - q_0 \right)^2}{\sum_K \lvert K \rvert \, q_0^2} \right)^{1/2},$$

and $e_v$ like $e_t$ on the vessel grid, where $\bar{q}_K$ is the computed exchange over the vessel control
volume $K$ divided by its length. The executable fails if the exchange that leaves the vessel and the exchange
that enters the tissue differ by more than $10^{-8} q_0 L$.

**Results**

The tables and figures are produced by

```bash
python3 run_benchmark.py
```

in the build directory `test/multidomain/embedded/1d3d/1p_1p_tip`. Passing `doc/doxygen/images` of the DuMux
source directory as argument regenerates the figures of this page.

![Errors of the coupling methods under grid refinement](embedded1d3d_tip_convergence.png)

| $h/R$ | $e_t$ LS | rate | $e_t$ CSS | rate | $e_t$ DS | rate | $e_t$ DS, box | rate |
|-------|----------|------|-----------|------|----------|------|---------------|------|
| 2 | $2.19 \times 10^{-2}$ | - | $1.90 \times 10^{-2}$ | - | $1.31 \times 10^{-2}$ | - | $4.84 \times 10^{-2}$ | - |
| 1 | $2.14 \times 10^{-2}$ | 0.04 | $3.01 \times 10^{-2}$ | -0.66 | $9.67 \times 10^{-3}$ | 0.44 | $1.85 \times 10^{-2}$ | 1.39 |
| 0.5 | $8.99 \times 10^{-3}$ | 1.25 | $1.17 \times 10^{-2}$ | 1.36 | $2.56 \times 10^{-3}$ | 1.92 | $4.55 \times 10^{-3}$ | 2.02 |
| 0.25 | $3.27 \times 10^{-3}$ | 1.46 | $4.42 \times 10^{-3}$ | 1.41 | $5.97 \times 10^{-4}$ | 2.10 | $1.15 \times 10^{-3}$ | 1.99 |

| $h/R$ | $e_v$ LS | rate | $e_v$ CSS | rate | $e_v$ DS | rate | $e_v$ DS, box | rate |
|-------|----------|------|-----------|------|----------|------|---------------|------|
| 2 | $9.06 \times 10^{-3}$ | - | $9.64 \times 10^{-3}$ | - | $2.14 \times 10^{-3}$ | - | $9.96 \times 10^{-3}$ | - |
| 1 | $4.70 \times 10^{-3}$ | 0.95 | $4.20 \times 10^{-3}$ | 1.20 | $1.70 \times 10^{-3}$ | 0.33 | $4.38 \times 10^{-3}$ | 1.18 |
| 0.5 | $2.05 \times 10^{-3}$ | 1.20 | $1.25 \times 10^{-3}$ | 1.75 | $4.71 \times 10^{-4}$ | 1.85 | $1.07 \times 10^{-3}$ | 2.03 |
| 0.25 | $7.48 \times 10^{-4}$ | 1.46 | $3.33 \times 10^{-4}$ | 1.90 | $1.07 \times 10^{-4}$ | 2.14 | $2.68 \times 10^{-4}$ | 2.00 |

| $h/R$ | $e_q$ LS | rate | $e_q$ CSS | rate | $e_q$ DS | rate | $e_q$ DS, box | rate |
|-------|----------|------|-----------|------|----------|------|---------------|------|
| 2 | $4.07 \times 10^{-2}$ | - | $4.33 \times 10^{-2}$ | - | $2.53 \times 10^{-2}$ | - | $1.15 \times 10^{-1}$ | - |
| 1 | $4.64 \times 10^{-2}$ | -0.19 | $4.37 \times 10^{-2}$ | -0.01 | $2.09 \times 10^{-2}$ | 0.27 | $4.50 \times 10^{-2}$ | 1.35 |
| 0.5 | $1.92 \times 10^{-2}$ | 1.28 | $1.29 \times 10^{-2}$ | 1.76 | $5.67 \times 10^{-3}$ | 1.88 | $1.12 \times 10^{-2}$ | 2.01 |
| 0.25 | $6.83 \times 10^{-3}$ | 1.49 | $3.51 \times 10^{-3}$ | 1.87 | $1.35 \times 10^{-3}$ | 2.07 | $2.79 \times 10^{-3}$ | 2.00 |

All methods converge towards their exact solutions; the closed tip does not lead to an error that stagnates
under refinement. The kernel method converges with approximately second order in all three quantities, with the
cell-centered and with the box method in the tissue. The line source method converges with an order between one
and two in all three quantities, the surface source method with an order between one and two in the tissue
pressure and with approximately second order in the vessel pressure and the exchange rate. These are the orders
of @ref benchmark-embedded-1d3d-vessel; the exchange error defined there uses the integrated exchange per control
volume, which adds one half to its observed orders.

![Tissue pressure in the plane through the vessel axis and vessel quantities along the vessel for h/R = 0.5](embedded1d3d_tip_solution.png)

The figure above shows the tissue pressure of the kernel method on the cells next to the plane $y = 0$ with
contour lines of the exact solution (dashed), the outline of the vessel (white), and the error of the vessel
pressure and the mean exchange rate per cell along the vessel for $h/R = 0.5$. The vessel pressure error is
smallest next to the inflow, where the pressure is prescribed, and largest at the tip. The mean exchange rate is
below $q_0$ in all cells. In the cells next to the tip it moves towards $q_0$ for the line source, the surface
source and the box variant of the kernel method, and away from $q_0$ for the cell-centered kernel method.

The test `test_md_embedded_1d3d_1p1p_tip_convergence` runs the script and fails if a rate between the two
finest grids is below 1.35 (LS), 1.3 for $e_t$ and 1.75 for $e_v$ and $e_q$ (CSS), or 1.9 (DS, both
schemes). The tests `test_md_embedded_1d3d_1p1p_tip_tpfatpfa_*` and `test_md_embedded_1d3d_1p1p_tip_boxtpfa_kernel`
run the setup of `params.input` with $c = 20$.
