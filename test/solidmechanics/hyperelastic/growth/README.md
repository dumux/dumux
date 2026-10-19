# Growth with Exact Solutions {#benchmark-hyperelastic-growth}

Growth, swelling and shrinkage enter the hyperelastic model through the multiplicative
decomposition \f$ \mathbf{F} = \mathbf{F}_\mathrm{e}\mathbf{F}_\mathrm{g} \f$ of the deformation
gradient into an elastic part and a growth tensor \f$ \mathbf{F}_\mathrm{g} \f$ @cite Rodriguez1994.
For a strain energy density \f$ \psi(\mathbf{F}_\mathrm{e}) \f$ per unit grown volume, the first
Piola-Kirchhoff stress is

$$\mathbf{P} = J_\mathrm{g}\, \mathbf{P}_\mathrm{e}(\mathbf{F}\mathbf{F}_\mathrm{g}^{-1})\,
\mathbf{F}_\mathrm{g}^{-T},
\quad J_\mathrm{g} = \det \mathbf{F}_\mathrm{g},$$

which the problems evaluate with `Growth::firstPiolaKirchhoffStressTensor` from
`dumux/solidmechanics/hyperelastic/growth.hh`. The benchmarks below have exact solutions: the
strip tests the kinematics of growth at large rotations, the two-layer cylinder and sphere the
stress pull-back with residual stress. The material is compressible neo-Hookean,
\f$ \mathbf{P}_\mathrm{e} = \mu \mathbf{F}_\mathrm{e} + (\lambda \ln J_\mathrm{e} - \mu)\mathbf{F}_\mathrm{e}^{-T} \f$,
with shear modulus \f$ \mu \f$ and Lamé parameter \f$ \lambda \f$. The growth is raised in steps
from zero, each solved with Newton's method from the previous solution. Each test runs three
nested meshes and checks the convergence rates of the \f$ L^2 \f$ and \f$ H^1 \f$ errors of the
displacement (`convergencetest.py`). The unit test `test_growth` compares the stress with
central differences of \f$ J_\mathrm{g}\psi(\mathbf{F}\mathbf{F}_\mathrm{g}^{-1}) \f$ for each
growth tensor and checks that \f$ \mathbf{F} = \mathbf{R}\mathbf{F}_\mathrm{g} \f$ with a
rotation \f$ \mathbf{R} \f$ is stress free.

`run_benchmark.py` runs all convergence studies of this page in the build directory
`test/solidmechanics/hyperelastic/growth` (about three minutes) and prints the tables; with the
directory `doc/doxygen/images` of the DuMux sources as argument, it also writes the figures.

## Strip bent by growth {#benchmark-hyperelastic-growth-strip}

The strip \f$ [0, L] \times [0, h] \f$ (3D: \f$ [0, L] \times [0, b] \times [0, h] \f$) is
clamped at \f$ X = 0 \f$ and traction free elsewhere. It grows along its axis with a stretch
that varies linearly with the distance \f$ T \f$ from the mid-surface,
\f$ \mathbf{F}_\mathrm{g} = \mathbf{I} + \kappa T\, \mathbf{e}_X \otimes \mathbf{e}_X \f$. This
growth is compatible: the deformation

$$x_X = (\kappa^{-1} + T)\sin(\kappa X), \quad
x_T = \tfrac{h}{2} + (\kappa^{-1} + T)\cos(\kappa X) - \kappa^{-1}$$

bends the strip into a circular arc of radius \f$ \kappa^{-1} \f$ and has
\f$ \mathbf{F} = \mathbf{R}\mathbf{F}_\mathrm{g} \f$, so it is stress free for any objective
material law and satisfies the boundary conditions. With \f$ L = 1 \f$, \f$ h = 0.1 \f$,
\f$ b = 0.1 \f$ and \f$ \kappa = \pi \f$ the strip closes to a half circle; the curvature is
raised in eight steps, \f$ \mu = 1 \f$, \f$ \lambda = 1.5 \f$.

![The strip of length L and thickness h, clamped at X = 0, grows along X with the stretch 1 + kappa T, which increases with the distance T from the mid-surface (dash-dotted). The exact solution bends it into a half circle of radius 1/kappa (shaded, dashed outline).](growth_strip_setup.svg){html: width=80%}

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

![Cross section of the two-layer cylinder or sphere: the outer layer (shaded) grows, the inner one does not; the layer interface is dashed. The quarter has symmetry conditions on the coordinate axes and traction-free inner and outer surfaces.](growth_shell_setup.svg){html: width=70%}

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

![Reference meshes (gray) and deformed meshes (black) of the strip and of the cylinder with the box method, and the exact boundaries and layer interface of the deformed bodies (dashed).](growth_shapes.png)

![L2 and H1 errors of the displacement against the largest element diameter, for the strip and for the cylinder and sphere, with all discretizations of this page including the mixed formulation below; each discretization has its own marker and line style, and the triangles show the slopes of the orders written next to them.](growth_convergence.png)

## Mixed displacement-pressure formulation {#benchmark-hyperelastic-growth-mixed}

The strip and the cylinder also run in the mixed formulation of @ref HyperelasticVolIso. The
neo-Hookean law is split as
\f$ \mathbf{P}_\mathrm{e} = \mu(\mathbf{F}_\mathrm{e} - \mathbf{F}_\mathrm{e}^{-T})
+ p_\mathrm{e} J_\mathrm{e}\mathbf{F}_\mathrm{e}^{-T} \f$ with
\f$ p_\mathrm{e} = \lambda \ln J_\mathrm{e}/J_\mathrm{e} \f$. The first part is pulled back with
`Growth::firstPiolaKirchhoffStressTensor`; the second becomes \f$ p J \mathbf{F}^{-T} \f$ with
\f$ p = \lambda \ln J_\mathrm{e}/J_\mathrm{e} \f$ and \f$ J_\mathrm{e} = J/J_\mathrm{g} \f$, which
the pressure problem returns from `volumetricPressure(J, elemDisc, ipGlobal)`. Both formulations
solve the same problem, so the exact solutions above apply. The tests use a nearly incompressible
material and report the \f$ L^2 \f$ error of the pressure besides the displacement errors.

In the strip, \f$ \lambda = 10^4 \f$ (\f$ \mu = 1 \f$). The exact solution is stress free, so
\f$ p = 0 \f$, while \f$ J = J_\mathrm{g} = 1 + \kappa T \f$ varies within the elements.

| Test | \f$ h_\mathrm{max} \f$ | \f$ L^2 \f$ error | rate | \f$ H^1 \f$ error | rate | pressure \f$ L^2 \f$ error | rate |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `test_hyperelastic_growth_bending_strip_mixed_p1b_2d` (MINI) | 0.0354 | 2.16e-4 | | 5.96e-3 | | 2.59e-3 | |
| | 0.0177 | 4.60e-5 | 2.23 | 2.91e-3 | 1.04 | 9.23e-4 | 1.49 |
| | 0.0088 | 1.00e-5 | 2.20 | 1.44e-3 | 1.02 | 3.26e-4 | 1.50 |
| `test_hyperelastic_growth_bending_strip_mixed_q2_2d` (Q2–Q1) | 0.0354 | 3.91e-5 | | 2.10e-4 | | 8.25e-5 | |
| | 0.0177 | 1.05e-5 | 1.90 | 5.55e-5 | 1.92 | 2.06e-5 | 2.00 |
| | 0.0088 | 2.66e-6 | 1.97 | 1.41e-5 | 1.98 | 5.16e-6 | 2.00 |

In the cylinder, \f$ \lambda = 1.5 \cdot 10^4 \f$ in the inner and \f$ 4.5 \cdot 10^4 \f$ in the
outer layer. The exact pressure jumps at the layer interface, which the continuous pressure
approximates with an \f$ L^2 \f$ error of order 1/2, and the displacement errors converge with
orders 1 and 1/2.

| Test | \f$ h_\mathrm{max} \f$ | \f$ L^2 \f$ error | rate | \f$ H^1 \f$ error | rate | pressure \f$ L^2 \f$ error | rate |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `test_hyperelastic_growth_layered_shell_mixed_p1b_2d` (MINI) | 0.1079 | 1.27e-3 | | 2.87e-2 | | 7.03e-2 | |
| | 0.0545 | 6.00e-4 | 1.10 | 1.98e-2 | 0.55 | 4.97e-2 | 0.51 |
| | 0.0274 | 2.90e-4 | 1.06 | 1.38e-2 | 0.52 | 3.51e-2 | 0.50 |
| `test_hyperelastic_growth_layered_shell_mixed_q2_2d` (Q2–Q1) | 0.1079 | 1.38e-3 | | 2.87e-2 | | 7.11e-2 | |
| | 0.0545 | 6.73e-4 | 1.06 | 2.02e-2 | 0.51 | 5.02e-2 | 0.51 |
| | 0.0274 | 3.32e-4 | 1.03 | 1.43e-2 | 0.50 | 3.55e-2 | 0.50 |

## Onset of growth-induced wrinkling {#benchmark-hyperelastic-growth-wrinkling}

A film of thickness \f$ h = 1 \f$ on a substrate of depth \f$ H = 30 \f$ grows along its length,
\f$ \mathbf{F}_\mathrm{g} = \mathbf{I} + (\vartheta - 1)\mathbf{e}_x \otimes \mathbf{e}_x \f$, in plane
strain. The substrate is clamped at the bottom, the sides \f$ x = 0 \f$ and \f$ x = L = 9 \f$ slide,
the top is traction free. Film and substrate are compressible neo-Hookean with \f$ \mu = 10 \f$ and
\f$ \mu = 1 \f$, \f$ \lambda = 1.5\mu \f$. The flat state, unstrained substrate and film compressed
along \f$ x \f$, solves the problem for every \f$ \vartheta \f$, and is exactly represented by the
discretization since its displacement is piecewise linear in \f$ y \f$ with a kink at a grid line.
It loses stability where the Jacobian becomes singular: the test raises \f$ \vartheta \f$ along the
flat solution until the smallest real eigenvalue of the Jacobian, computed with
`ShiftInvertEigenSolver`, changes sign, and locates the critical stretch by the Illinois method.

![A film of thickness h grows along x on a substrate of depth H, drawn shortened at the break; the substrate is clamped at the bottom and its sides slide, the top is traction free. The dashed curve above the film is the shape of the onset mode with two half wavelengths over the length L.](growth_wrinkling_setup.svg){html: width=80%}

The reference solves the problem linearized about the flat state (`reference.py`): in each layer a
perturbation \f$ u_x = U(y)\sin kx \f$, \f$ u_y = V(y)\cos kx \f$ satisfies a linear system of
ordinary differential equations with constant coefficients, the incremental moduli
\f$ \partial\mathbf{P}/\partial\mathbf{F} \f$ of the layer; propagating it from the clamped bottom
to the top gives the critical stretch as a root of a determinant. The sliding sides admit
\f$ k = m\pi/L \f$, and \f$ m = 2 \f$ wrinkles first, at \f$ \vartheta_c = 1.1430689262 \f$, before
\f$ m = 3 \f$ at 1.1634791. The propagator reproduces the stiffness
\f$ 4\mu k(1-\nu)/(3-4\nu) \f$ of a bonded half-space to \f$ 3 \cdot 10^{-11} \f$, and for a stiff
film on a deep substrate the critical strain approaches the classical asymptote with the bonded
substrate stiffness, with relative deviations of 3.2 %, 2.2 % and 1.4 % at stiffness ratios
1000, 3000 and 10000.

The meshes are tensor-product grids with the interface at a grid line, graded towards the film,
with 18 × (12 + 2) to 72 × (48 + 8) cells.

| Test | \f$ h_\mathrm{max} \f$ | \f$ \vartheta_c \f$ | error | rate | half wavelengths |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `test_hyperelastic_growth_wrinkling_box` | 4.84 | 1.1623740 | 1.93e-2 | | 2 |
| | 2.50 | 1.1478942 | 4.83e-3 | 2.10 | 2 |
| | 1.27 | 1.1442753 | 1.21e-3 | 2.05 | 2 |
| `test_hyperelastic_growth_wrinkling_pq2hybrid` (Q2) | 4.84 | 1.1448514 | 1.78e-3 | | 2 |
| | 2.50 | 1.1435616 | 4.93e-4 | 1.95 | 2 |
| | 1.27 | 1.1431998 | 1.31e-4 | 1.96 | 2 |

### Post-buckling branch {#benchmark-hyperelastic-growth-wrinkling-branch}

Beyond the onset, the wrinkled branch is traced by prescribing the vertical displacement of the top
corner at \f$ x = 0 \f$, a distance \f$ a \f$ above the flat state, and finding the growth stretch at
which its reaction vanishes by the secant method. Each amplitude starts from the previous state
moved along the eigenvector of the onset. The test `test_hyperelastic_growth_wrinkling_branch_box`
traces \f$ a = 0.01 \f$ to \f$ 0.16 \f$ on the coarsest mesh, fits
\f$ \vartheta = \vartheta_0 + c_2 a^2 + c_4 a^4 \f$ and checks that the intercept
\f$ \vartheta_0 \f$ agrees with the critical stretch from the eigenvalues to \f$ 10^{-6} \f$, that
\f$ c_2 > 0 \f$ (a supercritical pitchfork) and that \f$ \vartheta - \vartheta_c \f$ grows with
the exponent 2 between \f$ a = 0.01 \f$ and \f$ 0.02 \f$ to within 0.02. On the meshes of the
onset test, the intercept deviates from the critical stretch by \f$ 1.2 \cdot 10^{-7} \f$ to
\f$ 2.8 \cdot 10^{-7} \f$, and the exponent on the coarsest mesh is 2.0009.

| Discretization | \f$ c_2 \f$, coarsest mesh | middle mesh | finest mesh |
| :--- | ---: | ---: | ---: |
| box | 0.11843 | 0.12092 | 0.12140 |
| Q2 hybrid | 0.11854 | 0.12075 | 0.12133 |

`wrinkling/run_benchmark.py` runs the onset and the branch with the box and the Q2 hybrid scheme on
the meshes above in the build directory `test/solidmechanics/hyperelastic/growth/wrinkling` (about
two minutes) and prints the tables of this section; with the directory `doc/doxygen/images` of the
DuMux sources as argument, it also writes the figure.

![The film (shaded) and the top of the substrate on the box mesh with 36 cells along x at the largest traced amplitude, to scale; the growth stretch above the critical one on the branch against the amplitude for the box (hollow markers with lines) and the Q2 hybrid scheme (filled markers) on the meshes above, with the slope 2 of a supercritical pitchfork; the error of the critical stretch against the mesh size.](growth_wrinkling.png)
