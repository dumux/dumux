# Cook's Membrane, Linear Elasticity {#benchmark-cooks-membrane-elastic}

Cook's tapered panel @cite Cook1974 with a nearly incompressible linear elastic material, a standard
test for volumetric locking under combined bending and shear @cite Ong2015. The tests solve it with
the elastic model (@ref Elastic) assembled with local degrees of freedom, with three
control-volume finite-element schemes and two Galerkin finite-element schemes.

## Problem

The panel has the corners (0, 0), (48, 44), (48, 60) and (0, 44). The left edge \f$ x = 0 \f$ is
clamped, the right edge \f$ x = 48 \f$ carries the vertical traction \f$ t = 6.25 \f$ (a total load of
100), and the other edges are traction-free. Plane strain, Young's modulus \f$ E = 250 \f$ and
Poisson's ratio \f$ \nu = 0.4999 \f$, the setup of Ong et al. @cite Ong2015, who compare with the
results of Kasper and Taylor @cite Kasper2000. The quantity of interest is the vertical
displacement \f$ u_y \f$ of the upper right corner \f$ A = (48, 60) \f$. For the same load and
\f$ E \f$ with \f$ \nu = 0.4999999 \f$, Chen and Sukumar @cite Chen2024 quote the reference value
\f$ u_y(A) = 7.769 \f$.

## Discretization

| Test | Scheme |
| :--- | :--- |
| `test_cooks_membrane_box` | box (linear Lagrange basis, control volumes for all dofs) |
| `test_cooks_membrane_pq1bubble` | linear basis with element bubble, control volumes for all dofs |
| `test_cooks_membrane_pq2` | quadratic basis, control volumes for the vertex dofs and Galerkin test functions for the edge dofs |
| `test_cooks_membrane_pq1fe` | linear basis, Galerkin finite elements |
| `test_cooks_membrane_pq2fe` | quadratic basis, Galerkin finite elements |

The problem derives from `Experimental::ProblemWithSpatialParams`: the clamped edge is imposed as
constraints of the local degrees of freedom (`constraints()`) and the traction with
`boundaryFluxAtPos`. The mesh `cooksmembrane.msh` has 233 triangles; `Grid.Refinement` refines it
by conforming bisection, which roughly doubles the number of triangles per level.

## Results

Vertical displacement \f$ u_y(A) \f$ under refinement:

| Level | Box: dofs | Box | PQ2: dofs | PQ2 (control volumes) | PQ2 (finite elements) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 140 | 4.7086 | 512 | 7.5868 | 7.6491 |
| 1 | 303 | 7.0334 | 1140 | 7.5932 | 7.7000 |
| 2 | 632 | 7.4196 | 2429 | 7.7115 | 7.7306 |
| 3 | 1312 | 7.5569 | 5102 | 7.6931 | 7.7391 |
| 4 | 2709 | 7.6526 | 10,633 | 7.7376 | 7.7519 |
| 5 | 5577 | 7.6869 | 22,016 | 7.7332 | 7.7557 |
| 6 | 11,415 | 7.7204 | 45,249 | 7.7558 | 7.7619 |
| 7 | 23,237 | 7.7333 | 92,364 | 7.7530 | 7.7637 |
| 8 | 47,092 | 7.7477 | | | |

The PQ1-bubble and the PQ1 finite-element schemes give the same tip displacement as the box scheme on
every level computed (levels 0 to 4, agreement to \f$ 10^{-9} \f$). The values of the PQ2 control-volume
scheme alternate between odd and even levels; those of the PQ2 finite-element scheme increase
monotonically. Richardson extrapolation over every second level gives 7.759 (odd levels) and 7.766
(even levels) for the box scheme, 7.772 for the PQ2 control-volume scheme (odd levels) and 7.771 for
the PQ2 finite-element scheme (odd and even levels). These agree with the value 7.769 of Chen and
Sukumar @cite Chen2024 for \f$ \nu = 0.4999999 \f$ to 0.13 %. All schemes approach the limit slowly:
on the coarsest mesh the linear displacement elements lock (box: 4.71), and at the upper end of the
clamped edge, where clamped and free edge meet at 108°, the displacement behaves like \f$ r^{0.54} \f$
in the distance \f$ r \f$ from the corner (Williams' eigenvalue for \f$ \nu = 0.4999 \f$, plane strain).

The finite-difference Jacobian limits how close \f$ \nu \f$ can be to 0.5: at \f$ \nu = 0.4999999 \f$ the
Newton iteration of the PQ2 control-volume scheme diverges.

The CTests compare the displacement field on the base mesh with a reference.

```bash
./test_cooks_membrane_pq2 params.input -Problem.Name test_cooks_membrane_pq2 -Grid.Refinement 4
```
