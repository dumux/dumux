# Annular plate clamped at both edges

Uniform load on the annulus `a < r < b`, clamped at both edges, against the closed-form
axisymmetric deflection. The case has two boundary components that are entirely clamped,
so the gauge `phi = 0` cannot be prescribed on both: doing so drops the net flux of
`grad w - theta` through the hole, and the deflection converges to a wrong solution
(the test prints that error and the missing flux). The test then treats `phi` on the
inner edge as one unknown constant, determined by the vanishing summed compatibility
residual of the inner nodes, which is done by a secant in that constant because the
problem is linear in it, and checks second-order convergence of `w` and `phi` against the
closed form with `phi = -D (Delta w - Delta w(b))`. The constant is printed next to its
exact value `D (Delta w(b) - Delta w(a))`.

The third variant printed, the compatibility equation assembled naturally at every inner
node, converges in `w` but not in `phi`.

Run from the build directory:

```bash
ctest --output-on-failure -R '^test_kirchhoff_love_plate_annulus_clamped$'
```

The test checks three mesh sizes using `annulus.geo`. The VTK output contains the
solution with the unknown constant potential on the inner edge.
