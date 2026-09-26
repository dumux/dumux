# Kirchhoff-Love plate with free corners

This test measures the jump in twisting moment where two free edges meet. For a
convex corner without an applied point force, the Kirchhoff corner condition is
`R = M_ns(edge 1) - M_ns(edge 2) = 0`. With the tangent convention `s = J n`, the
twisting moments on the two edges of an axis-aligned right-angle corner are
`-M_xy` and `+M_xy`. The mixed formulation imposes `psi = M_ns` on each free edge;
continuity of `psi` supplies the corner condition.

The cases are a square, an L-shaped plate, and a right-triangular wedge. The square
and L-shape are clamped at `x = 0` and `y = 0`; the wedge is clamped only at `y = 0`.
All remaining edges are free. The test checks convex free-free corners, including
the acute corner of the wedge. It does not assess the singular re-entrant corner
of the L-shape. Clamped-free junctions are reported without a zero-force condition.

The moment tensor at each corner is averaged over its incident elements. The
reported ratio is `abs(R) / (D (1 - nu) max(abs(w)) / L^2)`. Across mesh sizes
`0.06`, `0.03`, and `0.018`, the finest ratio must be at most 60% of the coarsest
ratio at each tested free-free corner. Missing or nonfinite measurements fail the
test. The acceptance criterion measures reduction over these three meshes.

Run from the build directory:

```bash
ctest --output-on-failure -R '^test_kirchhoff_love_free_corner$'
```
