# Growth with Exact Solutions {#benchmark-hyperelastic-growth}

Growth, swelling and shrinkage enter the hyperelastic model through the multiplicative
decomposition \f$ \mathbf{F} = \mathbf{F}_\mathrm{e}\mathbf{F}_\mathrm{g} \f$ of the deformation
gradient into an elastic part and a growth tensor \f$ \mathbf{F}_\mathrm{g} \f$ @cite Rodriguez1994.
For a strain energy density \f$ \psi(\mathbf{F}_\mathrm{e}) \f$ per unit grown volume, the first
Piola-Kirchhoff stress is
\f[
    \mathbf{P} = J_\mathrm{g}\, \mathbf{P}_\mathrm{e}(\mathbf{F}\mathbf{F}_\mathrm{g}^{-1})\, \mathbf{F}_\mathrm{g}^{-T},
    \quad J_\mathrm{g} = \det \mathbf{F}_\mathrm{g},
\f]
which the problems evaluate with `Growth::firstPiolaKirchhoffStressTensor` from
`dumux/solidmechanics/hyperelastic/growth.hh`. The two benchmarks below have exact solutions:
the first tests the kinematics of growth at large rotations, the second the stress pull-back
with residual stress. The material is compressible neo-Hookean,
\f$ \mathbf{P}_\mathrm{e} = \mu \mathbf{F}_\mathrm{e} + (\lambda \ln J_\mathrm{e} - \mu)\mathbf{F}_\mathrm{e}^{-T} \f$,
with shear modulus \f$ \mu \f$ and Lamé parameter \f$ \lambda \f$. The growth is raised in steps
from zero, each solved with Newton's method from the previous solution. Each test runs three
nested meshes and checks the convergence rates of the \f$ L^2 \f$ and \f$ H^1 \f$ errors of the
displacement (`convergencetest.py`). The unit test `test_growth` compares the stress with
central differences of \f$ J_\mathrm{g}\psi(\mathbf{F}\mathbf{F}_\mathrm{g}^{-1}) \f$ for each
growth tensor and checks that \f$ \mathbf{F} = \mathbf{R}\mathbf{F}_\mathrm{g} \f$ with a
rotation \f$ \mathbf{R} \f$ is stress free.

## Strip bent by growth {#benchmark-hyperelastic-growth-strip}

The strip \f$ [0, L] \times [0, h] \f$ (3D: \f$ [0, L] \times [0, b] \times [0, h] \f$) is
clamped at \f$ X = 0 \f$ and traction free elsewhere. It grows along its axis with a stretch
that varies linearly with the distance \f$ T \f$ from the mid-surface,
\f$ \mathbf{F}_\mathrm{g} = \mathbf{I} + \kappa T\, \mathbf{e}_X \otimes \mathbf{e}_X \f$. This
growth is compatible: the deformation
\f[
    x_X = (\kappa^{-1} + T)\sin(\kappa X), \quad
    x_T = \tfrac{h}{2} + (\kappa^{-1} + T)\cos(\kappa X) - \kappa^{-1}
\f]
bends the strip into a circular arc of radius \f$ \kappa^{-1} \f$ and has
\f$ \mathbf{F} = \mathbf{R}\mathbf{F}_\mathrm{g} \f$, so it is stress free for any objective
material law and satisfies the boundary conditions. With \f$ L = 1 \f$, \f$ h = 0.1 \f$,
\f$ b = 0.1 \f$ and \f$ \kappa = \pi \f$ the strip closes to a half circle; the curvature is
raised in eight steps, \f$ \mu = 1 \f$, \f$ \lambda = 1.5 \f$.

| Test | Discretization | Mesh |
| :--- | :--- | :--- |
| `test_hyperelastic_growth_bending_strip_box_2d` | box | 40 × 4 to 160 × 16 squares |
| `test_hyperelastic_growth_bending_strip_pq2hybrid_2d` | Q2, hybrid control-volume / finite-element scheme | 40 × 4 to 160 × 16 squares |
| `test_hyperelastic_growth_bending_strip_pq2fe_2d` | Q2, finite elements | 40 × 4 to 160 × 16 squares |
| `test_hyperelastic_growth_bending_strip_box_3d` | box | 20 × 2 × 2 to 80 × 8 × 8 cubes |

| Test | \f$ h_\mathrm{max} \f$ | \f$ L^2 \f$ error | rate | \f$ H^1 \f$ error | rate |
| :--- | ---: | ---: | ---: | ---: | ---: |
| box, 2D | 0.0354 | 5.54e-3 | | 2.57e-2 | |
| | 0.0177 | 1.47e-3 | 1.91 | 7.49e-3 | 1.78 |
| | 0.0088 | 3.74e-4 | 1.98 | 2.46e-3 | 1.61 |
| Q2 hybrid, 2D | 0.0354 | 6.88e-5 | | 3.49e-4 | |
| | 0.0177 | 1.83e-5 | 1.91 | 9.27e-5 | 1.91 |
| | 0.0088 | 4.68e-6 | 1.97 | 2.37e-5 | 1.97 |
| Q2 finite elements, 2D | 0.0354 | 2.94e-6 | | 7.60e-5 | |
| | 0.0177 | 1.90e-7 | 3.95 | 1.85e-5 | 2.04 |
| | 0.0088 | 1.26e-8 | 3.91 | 4.59e-6 | 2.01 |
| box, 3D | 0.0866 | 6.33e-3 | | 2.81e-2 | |
| | 0.0433 | 1.93e-3 | 1.72 | 8.86e-3 | 1.67 |
| | 0.0217 | 5.10e-4 | 1.92 | 2.55e-3 | 1.80 |

## Growing two-layer cylinder and sphere {#benchmark-hyperelastic-growth-shell}

A hollow cylinder in plane strain (2D) or a hollow sphere (3D) consists of the layers
\f$ R_i \leq R < R_c \f$ and \f$ R_c < R \leq R_o \f$ with \f$ R_i = 0.5 \f$, \f$ R_c = 0.8 \f$,
\f$ R_o = 1 \f$. The outer layer grows isotropically,
\f$ \mathbf{F}_\mathrm{g} = 1.2\, \mathbf{I} \f$, the inner one does not; the moduli are
\f$ \mu = 1 \f$, \f$ \lambda = 1.5 \f$ inside and \f$ \mu = 3 \f$, \f$ \lambda = 4.5 \f$ outside.
Both surfaces are traction free. The growth mismatch leaves residual stress, and the solution
depends on the factor \f$ J_\mathrm{g}\mathbf{F}_\mathrm{g}^{-T} \f$ of the stress. A quarter
(2D) or an eighth (3D) is computed with symmetry conditions on a structured mesh whose layer
interface is a grid surface: radial lines over the quarter circle in 2D, the three faces of an
equiangular cubed sphere in 3D.

The exact solution is radially symmetric, \f$ \mathbf{x} = r(R)\, \mathbf{X}/R \f$. In a layer
with growth stretch \f$ \vartheta \f$ and \f$ n = d - 1 \f$ hoop directions, the principal stresses
are \f$ P_i = \vartheta^n (\mu a_i + (\lambda \ln J_\mathrm{e} - \mu)/a_i) \f$ with
\f$ a_r = r'/\vartheta \f$, \f$ a_\theta = r/(R\vartheta) \f$ and
\f$ J_\mathrm{e} = a_r a_\theta^n \f$, and radial equilibrium reads
\f$ P_r' + n(P_r - P_\theta)/R = 0 \f$. The test integrates this for \f$ (r, P_r) \f$, which is
continuous across the interface, with the classical Runge-Kutta method and finds
\f$ r(R_i) \f$ by shooting such that both surfaces are traction free (`radialreference.hh`).
An independent solution of the displacement form with a collocation boundary value solver
agrees to \f$ 10^{-11} \f$. The current radii of the inner surface, the interface and the outer
surface are 0.601108, 0.891269 and 1.138978 for the cylinder and 0.625936, 0.904834 and
1.157005 for the sphere.

| Test | Discretization | Mesh (radial inner × radial outer × angular cells) |
| :--- | :--- | :--- |
| `test_hyperelastic_growth_layered_shell_box_2d` | box | 6 × 4 × 16 to 24 × 16 × 64 quadrilaterals |
| `test_hyperelastic_growth_layered_shell_pq2hybrid_2d` | Q2, hybrid | 6 × 4 × 16 to 24 × 16 × 64 quadrilaterals |
| `test_hyperelastic_growth_layered_shell_box_3d` | box | 3 × 2 × (3 × 4²) to 12 × 8 × (3 × 16²) hexahedra |

| Test | \f$ h_\mathrm{max} \f$ | \f$ L^2 \f$ error | rate | \f$ H^1 \f$ error | rate |
| :--- | ---: | ---: | ---: | ---: | ---: |
| box, cylinder | 0.1079 | 5.22e-5 | | 4.03e-3 | |
| | 0.0545 | 1.30e-5 | 2.03 | 2.02e-3 | 1.02 |
| | 0.0274 | 3.26e-6 | 2.02 | 1.01e-3 | 1.01 |
| Q2 hybrid, cylinder | 0.1079 | 6.70e-5 | | 7.98e-4 | |
| | 0.0545 | 1.68e-5 | 2.02 | 2.69e-4 | 1.59 |
| | 0.0274 | 4.21e-6 | 2.01 | 9.28e-5 | 1.55 |
| box, sphere | 0.3017 | 3.60e-4 | | 2.43e-2 | |
| | 0.1592 | 9.14e-5 | 2.15 | 6.58e-3 | 2.04 |
| | 0.0818 | 2.29e-5 | 2.08 | 3.29e-3 | 1.04 |

The elements have straight edges, so the curved surfaces are approximated by polygons. The Q2
hybrid scheme reaches \f$ H^1 \f$ rates of 1.59 and 1.55 on the cylinder against 1.97 on the
strip, whose boundaries are straight.
