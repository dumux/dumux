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

The 2D mesh `cooksmembrane.msh` has 233 triangles, the 3D mesh `cooksmembrane3d.msh` 714 tetrahedra.
`Grid.Refinement` refines them by conforming bisection, which roughly doubles the number of elements
per level (535 and 1166 triangles on levels 1 and 2). The full load is applied in one Newton solve
with line search for each \f$ \lambda \f$.

## Results {#benchmark-cooks-membrane-hyperelastic-results}

### 2D, plane strain {#benchmark-cooks-membrane-hyperelastic-2d}

![Deformed membrane on the finest Taylor–Hood mesh (refinement 4), coloured by the displacement magnitude; dashed: undeformed membrane](cooks_membrane_hyperelastic_deformation.png)

Tip displacement \f$ u_y(A) \f$ in mm for Poisson's ratios 0.35 (\f$ \lambda = 432.099 \f$ MPa) and 0.4999
(\f$ \lambda = 926,000 \f$ MPa), and the number of degrees of freedom (displacement and pressure, including
the constrained ones):

| Refinement | MINI, dofs | MINI, \f$ \nu = 0.35 \f$ | MINI, \f$ \nu = 0.4999 \f$ | Taylor–Hood, dofs | Taylor–Hood, \f$ \nu = 0.35 \f$ | Taylor–Hood, \f$ \nu = 0.4999 \f$ |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1979 | 10.468 | 8.870 | 2583 | 10.567 | 9.001 |
| 2 | 4228 | 10.529 | 8.949 | 5490 | 10.584 | 9.024 |
| 3 | 8894 | 10.554 | 8.978 | 11,516 | 10.585 | 9.032 |
| 4 | 18,559 | 10.575 | 9.011 | 23,975 | 10.591 | 9.042 |
| 5 | 38,457 | 10.582 | 9.022 | | | |

![Tip displacement versus the number of degrees of freedom; reference values of Schröder et al. (2021) for Poisson's ratio 0.35](cooks_membrane_hyperelastic_convergence.png)

Schröder et al. report 10.5953 mm (p-FEM with plane-strain constraints, 18,213 dofs) and 10.5919 mm
(quadratic triangles T2, 115,680 dofs) (Table 3 of @cite Schroeder2021). At Poisson's ratio 0.35, Taylor–Hood
on refinement 4 (23,975 dofs) is 0.004 mm below the p-FEM value and 0.001 mm below the T2 value, and MINI on
refinement 5 is 0.013 mm below the p-FEM value. At Poisson's ratio 0.4999, the tip displacement of both
elements increases with every refinement; on the finest meshes the two differ by 0.020 mm, and Taylor–Hood
changes by 0.010 mm from refinement 3 to 4.

All five values of \f$ \lambda \f$ on the finest meshes (mm):

| \f$ \lambda \f$ [MPa] | MINI, refinement 5 | Taylor–Hood, refinement 4 |
| ---: | ---: | ---: |
| 432.099 | 10.582 | 10.591 |
| 750 | 10.093 | 10.105 |
| 9260 | 9.136 | 9.154 |
| 92,600 | 9.033 | 9.052 |
| 926,000 | 9.022 | 9.042 |

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
