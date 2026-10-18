# Vessel Crossing a Tissue Block (1D-3D Embedded Coupling) {#benchmark-embedded-1d3d-vessel}

## Convergence of line source, surface source and kernel coupling towards exact solutions

**Problem description**

A straight vessel of radius $R$ with center-line $\Lambda = \{ (0, 0, s) : s \in [0, 1] \}$ crosses the
tissue block $\Omega = [-1, 1] \times [-1, 1] \times [0, 1]$. Fluid flows along the vessel and leaks into the
tissue. The vessel is described by a one-dimensional model on $\Lambda$ and the tissue by a three-dimensional
model on $\Omega$, which also covers the space occupied by the vessel (mixed-dimension embedded model),

$$-\frac{\mathrm{d}}{\mathrm{d}s} \left( \kappa_v \frac{\mathrm{d} p_v}{\mathrm{d}s} \right) = -q \quad \text{on } \Lambda,
\qquad
-\Delta p_t = q \, \delta_\Lambda \quad \text{in } \Omega,$$

where $p_v$ is the vessel pressure, $\kappa_v$ the axial conductivity of the vessel, $p_t$ the tissue pressure,
and $q$ the exchange rate per unit vessel length, positive from the vessel into the tissue. The permeability of
the tissue and the fluid viscosity are one. This is the single-vessel benchmark of Koch et al. (2020)
@cite Koch2020 (Section 4.2), there posed on the unit cube with the vessel at $x = y = 0.5$.

The distribution $\delta_\Lambda$ places the exchange in the tissue. The coupling methods differ in this
distribution and in the tissue pressure that enters the exchange rate. With $r = \sqrt{x^2 + y^2}$ the distance
from the center-line and the coefficient $\beta = 2\pi / (2\pi + \ln R)$,

| Method | Coupling mode | Source in the tissue | Exchange rate $q$ |
|--------|---------------|----------------------|-------------------|
| line source (LS) | `Embedded1d3dCouplingMode::Average` | on the center-line | $\beta (p_v - \bar{p}_t)$ |
| cylinder surface source (CSS) | `Embedded1d3dCouplingMode::Surface` | uniform on the vessel surface $r = R$ | $\beta (p_v - \bar{p}_t)$ |
| distributed source (DS) | `Embedded1d3dCouplingMode::Kernel` | uniform in the cylinder $r \leq \varrho$ | $\beta \xi (p_v - p_{t,0})$ |

where $\bar{p}_t$ is the average of $p_t$ over the circle of radius $R$ around the center-line, $p_{t,0}$ the
tissue pressure on the center-line, $\varrho$ the kernel radius and $\xi$ the flux scaling factor of the
kernel method. In the continuous model, $\xi = \left( 1 + \frac{\beta}{2\pi} \left( \ln \frac{\varrho}{R} - \frac12 \right) \right)^{-1}$
@cite Koch2020 (Eq. 2.12), so that $\beta \xi (p_v - p_{t,0}) = \beta (p_v - \bar{p}_t)$. The discrete method
evaluates $\xi$ at the mean distance between the vessel segment and the tissue cell (Eqs. 3.18 and 3.20).

![Setup in the plane y = 0 through the vessel axis with the boundary conditions (solid: prescribed tissue pressure, dashed: prescribed normal flux). Right: the source of each method in a cross-section of the vessel; the line source method averages the tissue pressure over the dashed circle of radius R.](embedded1d3d_vessel_setup.svg){html: width=100%}

**Exact solution**

With the vessel conductivity $\kappa_v = 1 + s + s^2/2$, the pressures

$$p_v = 1 + s,
\qquad
p_t^\mathrm{LS} = -\frac{1 + z}{2\pi} \ln r,
\qquad
p_t^\mathrm{CSS} = -\frac{1 + z}{2\pi} \ln \max(r, R),$$

$$p_t^\mathrm{DS} = \begin{cases}
-\frac{1 + z}{2\pi} \left( \frac{r^2}{2\varrho^2} + \ln \varrho - \frac12 \right) & r \leq \varrho, \\
-\frac{1 + z}{2\pi} \ln r & r > \varrho,
\end{cases}$$

solve the problem with the exchange rate $q = 1 + s$ for each method. The tissue pressures of the three methods
agree for $r \geq \max(R, \varrho)$. Because $q$ is linear in $s$, the tissue pressure in each plane
$z = \mathrm{const}$ is the two-dimensional potential of a point source.

| Boundary | Condition |
|----------|-----------|
| vessel, $s = 0$ and $s = 1$ | $p_v = 1$ and $p_v = 2$ |
| tissue, $x = \pm 1$ and $y = \pm 1$ | exact tissue pressure |
| tissue, $z = 0$ and $z = 1$ | normal flux of the exact tissue pressure |

**Parameters**

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Vessel radius (convergence study) | $R$ | 0.1 |
| Vessel radius (`params.input`) | $R$ | 0.3 |
| Kernel radius | $\varrho$ | $R$ (`MixedDimension.KernelWidthFactor = 1`) |
| Tissue permeability, fluid viscosity, fluid density | | 1 |
| Points of the circle average | | 256 (`MixedDimension.NumCircleSegments`) |
| Quadrature order on the vessel segments | | 2 (`MixedDimension.IntegrationOrder`) |

**Discretization**

The tissue grid has $c \times c \times c/2$ cells (`Dune::YaspGrid`), the vessel grid $c/2$ cells
(`Dune::FoamGrid`), so that both have the cell size $h = 2/c$. Both domains use the cell-centered
finite volume method with two-point flux approximation (@ref CCTpfaDiscretization) or the box method
(@ref BoxDiscretization). The coupling is implemented by the coupling managers of @ref EmbeddedCoupling.
The problem is linear and is solved with a single Newton step.

The test computes the errors

$$e_t = \frac{\left( \sum_K \lvert K \rvert \left( p_{t,K} - p_t(\mathbf{x}_K) \right)^2 \right)^{1/2}}{\sum_K \lvert K \rvert},
\qquad
e_v = \frac{\left( \sum_K \lvert K \rvert \left( p_{v,K} - p_v(\mathbf{x}_K) \right)^2 \right)^{1/2}}{\sum_K \lvert K \rvert},
\qquad
e_q = \frac{\left( \sum_K \left( Q_K - Q_{K,e} \right)^2 \right)^{1/2}}{\lvert \Lambda \rvert},$$

where the sums run over the cells $K$ of the tissue or the vessel grid with volume or length $\lvert K \rvert$ and
center $\mathbf{x}_K$, $p_{t,K}$ and $p_{v,K}$ are the discrete solutions evaluated at $\mathbf{x}_K$,
$Q_K$ is the computed exchange over the vessel control volume $K$, and $Q_{K,e} = \int_K (1 + s) \, \mathrm{d}s$.

**Results**

The convergence study runs the cell-centered schemes in both domains for $c = 10, 20, 40, 80$ and $R = 0.1$,
i.e. $h/R = 2, 1, 0.5, 0.25$:

```bash
python3 convergence.py
```

in the build directory `test/multidomain/embedded/1d3d/1p_1p`. Passing `doc/doxygen/images` of the DuMux
source directory as argument regenerates the figure of this page.

![Errors of the line source, surface source and kernel methods under grid refinement](embedded1d3d_vessel_convergence.png)

| $h/R$ | $e_t$ LS | rate | $e_t$ CSS | rate | $e_t$ DS | rate |
|-------|----------|------|-----------|------|----------|------|
| 2 | $4.06 \times 10^{-3}$ | - | $4.05 \times 10^{-3}$ | - | $1.32 \times 10^{-3}$ | - |
| 1 | $3.95 \times 10^{-3}$ | 0.04 | $4.33 \times 10^{-3}$ | -0.10 | $1.47 \times 10^{-3}$ | -0.15 |
| 0.5 | $1.68 \times 10^{-3}$ | 1.24 | $1.58 \times 10^{-3}$ | 1.46 | $3.77 \times 10^{-4}$ | 1.96 |
| 0.25 | $6.10 \times 10^{-4}$ | 1.46 | $5.65 \times 10^{-4}$ | 1.48 | $8.45 \times 10^{-5}$ | 2.16 |

| $h/R$ | $e_v$ LS | rate | $e_v$ CSS | rate | $e_v$ DS | rate |
|-------|----------|------|-----------|------|----------|------|
| 2 | $4.99 \times 10^{-3}$ | - | $4.98 \times 10^{-3}$ | - | $7.07 \times 10^{-4}$ | - |
| 1 | $3.91 \times 10^{-3}$ | 0.35 | $3.92 \times 10^{-3}$ | 0.34 | $1.48 \times 10^{-3}$ | -1.06 |
| 0.5 | $1.65 \times 10^{-3}$ | 1.24 | $1.18 \times 10^{-3}$ | 1.73 | $3.70 \times 10^{-4}$ | 1.99 |
| 0.25 | $5.98 \times 10^{-4}$ | 1.47 | $3.32 \times 10^{-4}$ | 1.83 | $8.11 \times 10^{-5}$ | 2.19 |

| $h/R$ | $e_q$ LS | rate | $e_q$ CSS | rate | $e_q$ DS | rate |
|-------|----------|------|-----------|------|----------|------|
| 2 | $3.98 \times 10^{-2}$ | - | $4.00 \times 10^{-2}$ | - | $5.64 \times 10^{-3}$ | - |
| 1 | $2.30 \times 10^{-2}$ | 0.79 | $2.29 \times 10^{-2}$ | 0.80 | $8.61 \times 10^{-3}$ | -0.61 |
| 0.5 | $6.96 \times 10^{-3}$ | 1.73 | $4.92 \times 10^{-3}$ | 2.22 | $1.55 \times 10^{-3}$ | 2.48 |
| 0.25 | $1.79 \times 10^{-3}$ | 1.96 | $9.82 \times 10^{-4}$ | 2.32 | $2.40 \times 10^{-4}$ | 2.69 |

Only the kernel method converges with approximately second order in the tissue pressure, and it has the smallest errors in
all three quantities on all grids, as reported by Koch et al. (2020). Between $h/R = 2$ and $h/R = 1$ the errors
hardly decrease; for $h > R$, the sources of all three methods lie in the cells adjacent to the center-line.
On uniform vessel grids, $Q_K - Q_{K,e} = h (\bar{q}_K - \bar{q}_{K,e})$ with the mean exchange rates
$\bar{q}_K$ per cell, so that $e_q = h^{1/2} \left( \sum_K h (\bar{q}_K - \bar{q}_{K,e})^2 \right)^{1/2} / \lvert \Lambda \rvert$
and the observed orders of $e_q$ are higher by one half than those of the discrete $L^2$ error of the mean
exchange rate.

The test `test_md_embedded_1d3d_1p1p_tpfatpfa_convergence` runs the script and fails if the mean rates of the
two finest refinements, $(e_t, e_v, e_q)$ per method, differ from the rates stored in the script by more than
0.1 in the Euclidean norm.

The grids of the convergence study have the center-line on cell faces and an integer ratio $R/h$ for $h \leq R$.
Without this alignment the errors do not decrease monotonically: with $R = 0.1$ and $c = 32, 64$
($h/R = 0.625, 0.3125$), $e_t$ changes from $3.49 \times 10^{-4}$ to $3.14 \times 10^{-4}$ (LS), from
$6.96 \times 10^{-4}$ to $6.69 \times 10^{-4}$ (CSS), and from $2.66 \times 10^{-4}$ to $1.39 \times 10^{-4}$ (DS).
With $R = 0.125$ on the same grids ($h/R = 0.5, 0.25$), the errors decrease again with the orders of the tables
above. Finer kernel integration (`MixedDimension.KernelIntegrationCRL` 0.003 instead
of 0.01) or more circle points (2048 instead of 256) change the errors at $c = 64$ by less than 3 %.

**Further tests**

- `test_md_embedded_1d3d_1p1p_*_average` compare the solution with reference files for all combinations of the
  cell-centered and the box method. The box tissue grid has an odd number of cells per direction, so that no
  vertex lies on the center-line, where $p_t^\mathrm{LS}$ is singular.
- `test_md_embedded_1d3d_1p1p_*_kernel` check that the kernel weights integrate to the vessel length.
- `test_md_embedded_1d3d_1p1p_tpfatpfa_line` runs `Embedded1d3dCouplingMode::Line`, which evaluates the tissue
  pressure on the center-line. The exact line source solution is singular there, so this mode has no exact
  solution in this setup.
- `test_md_embedded_1d3d_1p1p_boxtpfa_kernel_unstructured` runs the kernel method with the box method on
  tetrahedral meshes generated from `box.geo` with Gmsh, and fails if the rate of $e_t$ between the two
  finest meshes is below 1.5.

The benchmark @ref benchmark-embedded-1d3d-vessel-tip uses the same model with a vessel that ends inside the
tissue.
