# Cook's Membrane, Hyperelastic (SPP 1748) {#benchmark-cooks-membrane-hyperelastic}

Cook's tapered cantilever with a hyperelastic material, benchmark 1 (Sect. 2.3) of the selection of
solid mechanics benchmarks of the DFG priority programme SPP 1748, Schröder et al. (2021)
@cite Schroeder2021. The membrane combines bending and shear. Its compressibility is reduced in five
steps from Poisson's ratio 0.35 to 0.4999, which shows the volumetric locking of displacement-only
elements of low order. The tests solve it with the mixed displacement–pressure formulation with the
volumetric–isochoric split (@ref HyperelasticVolIso), in 2D (this folder) and in 3D
(`test/solidmechanics/hyperelastic_voliso/cooks_membrane_3d`).

## Problem {#benchmark-cooks-membrane-hyperelastic-problem}

![Geometry and boundary conditions of Cook's membrane, dimensions in mm](cooks_membrane_setup.svg){html: width=50%}

The membrane has the corners (0, 0), (48, 44), (48, 60) and (0, 44) mm and the thickness 1 mm. The
left edge \f$ x = 0 \f$ is clamped and the right edge \f$ x = 48 \f$ mm carries the shear traction
\f$ p_0 = 20 \f$ MPa in \f$ +y \f$ direction (dead load); the other edges are free. The material is
\f[
\psi_1 = \frac{\mu}{2}(\operatorname{tr}\mathbf{C} - 3) + \frac{\lambda}{4}(J^2 - 1) - \left(\frac{\lambda}{2} + \mu\right) \ln J,
\qquad
\mathbf{P} = \mu (\mathbf{F} - \mathbf{F}^{-\top}) + \frac{\lambda}{2}(J^2 - 1)\, \mathbf{F}^{-\top},
\f]
with \f$ \mu = 185.185 \f$ MPa and \f$ \lambda \f$ = 432.099, 750, 9260, 92,600 and 926,000 MPa (Poisson's ratio
0.35, 0.401, 0.4902, 0.499 and 0.4999); every test runs all five values. The quantity of interest is
the vertical displacement \f$ u_y \f$ of the upper right corner \f$ A = (48, 60) \f$ mm.

- 2D: plane strain.
- 3D: the front and back faces are free. Half of the thickness, \f$ 0 \leq z \leq 0.5 \f$ mm, is modelled
  with \f$ u_z = 0 \f$ on the mid-plane \f$ z = 0 \f$; \f$ A = (48, 60, 0) \f$ mm.

The plane-strain and the 3D values differ by about 7 % because the 3D membrane is not constrained
in thickness direction.

## Discretization {#benchmark-cooks-membrane-hyperelastic-discretization}

| Test | Displacement | Pressure |
| :--- | :--- | :--- |
| `test_hyperelastic_cooks_membrane_mini`, `test_hyperelastic_cooks_membrane_mini_3d` | P1 with bubble | P1 |
| `test_hyperelastic_cooks_membrane_taylorhood`, `test_hyperelastic_cooks_membrane_taylorhood_3d` | P2 | P1 |

The 2D mesh `test/solidmechanics/cooksmembrane.msh` has 233 triangles, the 3D mesh `cooksmembrane3d.msh` 714 tetrahedra.
`Grid.Refinement` refines them by conforming bisection, which roughly doubles the number of elements
per level (535 and 1166 triangles on levels 1 and 2). The full load is applied in one Newton solve
with line search for each \f$ \lambda \f$.

## Results {#benchmark-cooks-membrane-hyperelastic-results}

### 2D, plane strain {#benchmark-cooks-membrane-hyperelastic-2d}

![Deformed membrane on Taylor–Hood refinement 4 for Poisson's ratios 0.35 and 0.4999, coloured by the displacement magnitude; dashed: undeformed membrane](cooks_membrane_hyperelastic_deformation.png)

Tip displacement \f$ u_y(A) \f$ in mm for Poisson's ratios 0.35 (\f$ \lambda = 432.099 \f$ MPa) and 0.4999
(\f$ \lambda = 926,000 \f$ MPa), and the number of degrees of freedom (displacement and pressure, including
the constrained ones). The last row is the Richardson extrapolation of the last three even and the last
three odd refinement levels:

| Refinement | MINI, dofs | MINI, \f$ \nu = 0.35 \f$ | MINI, \f$ \nu = 0.4999 \f$ | Taylor–Hood, dofs | Taylor–Hood, \f$ \nu = 0.35 \f$ | Taylor–Hood, \f$ \nu = 0.4999 \f$ |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1979 | 10.4676 | 8.8703 | 2583 | 10.5672 | 9.0009 |
| 2 | 4228 | 10.5291 | 8.9490 | 5490 | 10.5837 | 9.0239 |
| 3 | 8894 | 10.5535 | 8.9784 | 11,516 | 10.5847 | 9.0315 |
| 4 | 18,559 | 10.5749 | 9.0106 | 23,975 | 10.5912 | 9.0415 |
| 5 | 38,457 | 10.5817 | 9.0222 | 49,609 | 10.5910 | 9.0455 |
| 6 | 79,085 | 10.5888 | 9.0363 | 101,913 | 10.5937 | 9.0500 |
| 7 | 161,493 | 10.5906 | 9.0409 | 207,965 | 10.5933 | 9.0525 |
| 8 | 327,992 | 10.5930 | 9.0475 | | | |
| extrapolated | | 10.5947–10.5948 | 9.0549–9.0563 | | 10.5947–10.5949 | 9.0577–9.0592 |

![Relative error of the tip displacement at Poisson's ratio 0.35 with respect to the p-FEM value of Schröder et al. (2021); dotted: deviation of the extrapolated value](cooks_membrane_hyperelastic_convergence.png)

Schröder et al. report 10.5953 mm (p-FEM with plane-strain constraints, 18,213 dofs) and 10.5919 mm
(quadratic triangles T2, 115,680 dofs) (Table 3 of @cite Schroeder2021). `Grid.Refinement` bisects every
triangle once per level, so the mesh size halves every second level and the tip displacement changes in
steps of two levels. At Poisson's ratio 0.35, the change over two levels decreases by a factor of 2.7 to 3.3
from one pair of levels to the next for both elements, roughly \f$ N^{-3/4} \f$ in the number of degrees of
freedom \f$ N \f$. Both elements extrapolate to 10.5947 to 10.5949 mm, 0.0004 to 0.0006 mm below the p-FEM
value and 0.0028 to 0.0030 mm above the T2 value. Taylor–Hood on refinement 6 (101,913 dofs) is 0.0016 mm
below the p-FEM value, T2 with 115,680 dofs 0.0034 mm.

At Poisson's ratio 0.4999, the factor is 2.0 to 2.3 and the values extrapolate to 9.0549 to 9.0592 mm.
The finest meshes are 0.06 % to 0.07 % (Taylor–Hood, refinement 7) and 0.08 % to 0.10 % (MINI, refinement 8)
below the extrapolated values of the same element.

All five values of \f$ \lambda \f$ on the finest meshes and extrapolated (mm):

| \f$ \lambda \f$ [MPa] | MINI, refinement 8 | Taylor–Hood, refinement 7 | extrapolated, both elements |
| ---: | ---: | ---: | ---: |
| 432.099 | 10.5930 | 10.5933 | 10.5947–10.5949 |
| 750 | 10.1082 | 10.1094 | 10.1109–10.1122 |
| 9260 | 9.1599 | 9.1642 | 9.1665–9.1705 |
| 92,600 | 9.0580 | 9.0629 | 9.0653–9.0696 |
| 926,000 | 9.0475 | 9.0525 | 9.0549–9.0592 |

Schröder et al. give no plane-strain reference for \f$ \lambda > 432.099 \f$ MPa.

### 3D {#benchmark-cooks-membrane-hyperelastic-3d}

Tip displacement \f$ u_y(A) \f$ in mm on the finest meshes, with the p-FEM reference values of Schröder
et al. (Table 9 of @cite Schroeder2021):

| \f$ \lambda \f$ [MPa] | reference | MINI, refinement 3 | Taylor–Hood, refinement 2 |
| ---: | ---: | ---: | ---: |
| 432.099 | 11.381 | 11.183 | 11.349 |
| 750 | 11.156 | 10.933 | 11.119 |
| 9260 | 10.783 | 10.420 | 10.735 |
| 92,600 | 10.746 | 10.332 | 10.698 |
| 926,000 | 10.741 | 10.321 | 10.694 |

At Poisson's ratio 0.4999, MINI reaches 95 %, 95 % and 96 % of the reference on refinement levels 1,
2 and 3, and Taylor–Hood 99.4 % and 99.6 % on levels 1 and 2.

## Running the tests {#benchmark-cooks-membrane-hyperelastic-running}

```bash
./test_hyperelastic_cooks_membrane_taylorhood params_mixed.input -Grid.Refinement 3
./test_hyperelastic_cooks_membrane_mini_3d params_mixed3d.input -Grid.Refinement 3
```

The 2D results and figures are reproduced in the build directory of the 2D tests with
(requires numpy, matplotlib and pyvista)

```bash
python3 plot_results.py ./test_hyperelastic_cooks_membrane_mini ./test_hyperelastic_cooks_membrane_taylorhood
```
