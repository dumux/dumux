# Incompressible Block (SPP 1748) {#benchmark-incompressible-block}

A nearly incompressible hyperelastic block under a partial surface load, benchmark 2 (Sect. 3) of
the selection of solid mechanics benchmarks of the DFG priority programme SPP 1748, Schröder et al.
(2021) @cite Schroeder2021. Two line singularities of the load intersect in a point, and the material
is close to incompressible, so displacement-only elements of low order lock. The tests solve it with
the mixed displacement–pressure formulation with the volumetric–isochoric split (@ref HyperelasticVolIso).

## Problem {#benchmark-incompressible-block-problem}

![Modelled quarter of the block with its boundary conditions, dimensions in mm](incompressible_block_setup.svg){html: width=50%}

The block has the width and length 100 mm and the height 50 mm. Because of the symmetry of the
loading, one quarter, \f$ [0, 50]^3 \f$ mm, is modelled:
- bottom \f$ z = 0 \f$: \f$ u_z = 0 \f$,
- top \f$ z = 50 \f$ mm: \f$ u_x = u_y = 0 \f$,
- symmetry planes \f$ x = 50 \f$ mm and \f$ y = 50 \f$ mm: \f$ u_x = 0 \f$ and \f$ u_y = 0 \f$, respectively,
- the part \f$ 25 \leq x, y \leq 50 \f$ mm of the top face (the quarter of the 50 mm × 50 mm loaded
  area at the centre of the block) carries the dead load \f$ q = 3 \f$ MPa in \f$ -z \f$ direction.

The strain energy density is
\f$ W = \frac{\mu}{2}(\operatorname{tr}\mathbf{C} - 3 - 2 \ln J) + \frac{\lambda}{2}(J - 1)^2 \f$
with \f$ \lambda = 499.92568 \f$ MPa and \f$ \mu = 1.61148 \f$ MPa (Poisson's ratio 0.4983). The quantity
of interest is the vertical displacement \f$ |u_z(P)| \f$ of the centre \f$ P = (50, 50, 50) \f$ mm of the
loaded area. Schröder et al. report a converged value of 20.0 mm, and for the trilinear
displacement element (H1) 7.78 mm (260 dofs), 13.17 mm (1800 dofs), 17.54 mm (13,328 dofs) and
19.52 mm (102,432 dofs) (Table 17 of @cite Schroeder2021).

## Discretization {#benchmark-incompressible-block-discretization}

| Test | Displacement | Pressure | Mesh |
| :--- | :--- | :--- | :--- |
| `test_hyperelastic_incompressible_block_p1b_p1` | P1 with bubble (MINI) | P1, box | tetrahedra, `block_delaunay.msh` |
| `test_hyperelastic_incompressible_block_p2_p1` | P2, hybrid control-volume / finite-element scheme | P1, box | tetrahedra, `block_delaunay.msh` |
| `test_hyperelastic_incompressible_block_p2fe_p1fe` | P2, finite elements | P1, finite elements | tetrahedra, `block_delaunay.msh` |
| `test_hyperelastic_incompressible_block_q2_q1` | Q2, hybrid | Q1, box | 4 × 4 × 4 hexahedra |

The load is applied in five steps (`LoadStepping.NSteps`).

## Results {#benchmark-incompressible-block-results}

Vertical displacement \f$ |u_z(P)| \f$ in mm:

| Test | dof positions (displacement + pressure) | \f$ \lvert u_z(P) \rvert \f$ |
| :--- | ---: | ---: |
| P1-bubble–P1 (MINI) | 818 | 19.26 |
| P2–P1, hybrid | 1112 | 20.02 |
| P2–P1, finite elements | 1112 | 19.97 |
| Q2–Q1, hybrid | 854 | 20.10 |

With about a thousand dof positions, the mixed formulations are within 3.7 % of the converged value
20.0 mm. The trilinear displacement element of Schröder et al. gives 13.17 mm with 1800 dofs.

![Deformed block after the last load step, Q2–Q1 on 8 × 8 × 8 hexahedra per quarter, coloured by the displacement magnitude; grey lines: undeformed block](incompressible_block_deformation.png)

On the 8 × 8 × 8 hexahedral mesh (5642 dof positions), Q2–Q1 gives \f$ |u_z(P)| \f$ = 20.14 mm.

The CTests compare the displacement field with a reference. The figure is produced in the build
directory with (requires numpy, matplotlib and pyvista)

```bash
python3 plot_results.py ./test_hyperelastic_incompressible_block_q2_q1
```
