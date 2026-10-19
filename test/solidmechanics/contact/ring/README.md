# Collapse of a Ring under External Pressure {#benchmark-contact-ring}

## Self-contact of a thin elastic ring buckled by pressure, compared with the elastica

**Problem description**

A thin elastic ring of midline radius $R$ and wall thickness $t$ is loaded by an external
pressure $p$ (plane strain). At $p = 3D/R^3$ it buckles into an elliptic shape, with the bending
stiffness $D = E t^3 / (12 (1 - \nu^2))$, Young's modulus $E$ and Poisson ratio $\nu$. Under
growing pressure the opposite walls approach each other, touch at a point and then along a
growing line. The benchmark computes this collapse with self-contact of the inner surface
(`Contact::PenaltyContact`, `Contact::ContactAssembler`) and compares the contact force and the
onset of line contact with the inextensible elastica @cite Flaherty1972.

The material is neo-Hookean, $\mathbf{P} = \mu (\mathbf{F} - \mathbf{F}^{-T}) + \lambda \ln J \,
\mathbf{F}^{-T}$, with the deformation gradient $\mathbf{F}$, $J = \det \mathbf{F}$ and the Lamé
parameters $\mu$, $\lambda$ of $E$ and $\nu$. The pressure follows the deformed outer surface: the
traction per reference area is $-p (1 - \epsilon \cos 2\varphi) \operatorname{cof}(\mathbf{F})
\mathbf{N}$ with the outer unit normal $\mathbf{N}$ in the reference configuration and the polar
angle $\varphi$. The small imperfection $\epsilon$ selects the elliptic mode with the minor axis
along $y$. The nodes on the coordinate axes stay on them, which removes the rigid motions without
restricting the deformation, which is symmetric about both axes. The inner surface is the contact
surface; the penalty law with stiffness $k$ and offset $\delta$ is described on the page
@ref benchmark-contact-hertz.

![Left: the ring of midline radius R and wall thickness t under the external pressure p (1 - epsilon cos 2 phi); the dash-dotted lines are the coordinate axes, on which the nodes stay. Right: the collapsed ring with its inner walls in contact along a flat segment of length 2b.](ring_setup.svg)

**Reference solution**

The ring is modelled as an inextensible elastica with bending stiffness $D$. The elliptic mode is
symmetric about both axes, so a quarter of the midline suffices, from the major axis to the minor
axis. With the tangent angle $\theta$ and the arc length $s$, moment balance gives

$$D \theta' = C + \frac{p}{2} \left( x^2 + y^2 \right) - \frac{Q}{2} x,$$

with the moment constant $C$ and the contact force $Q$ per unit length between the two halves of
the ring. Contact of the walls means that the midline reaches $y = h = t/2$ on the minor axis.

- Before contact, $Q = 0$, and the free elliptic branch starts at $p = 3 D/R^3$.
- At first contact the midline touches $y = h$ at a point, and $Q$ grows with $p$ as long as the
  curvature at the contact point stays negative.
- Line contact starts when this curvature vanishes. Then the midline is flat on $|x| \leq b$,
  with the pressure $p$ on the flat part and point forces at its ends.

`elastica.py` solves each branch by shooting with an adaptive Runge-Kutta method. For $h = 0$
it gives first contact at $p R^3/D = 5.246887$ and line contact at $10.336033$, which agree with
the values 5.247 and 10.34 of @cite Flaherty1972; @cite Djondjorov2011 confirm 5.247
independently. This benchmark uses $h/R = 0.02$.

**Parameters**

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| midline radius | $R$ | 1 | m |
| wall thickness | $t$ | 0.04 | m |
| Young's modulus | $E$ | 1 | Pa |
| Poisson ratio | $\nu$ | 0.3 | - |
| imperfection | $\epsilon$ | 0.01 | - |
| penalty stiffness | $k$ | 10 | Pa/m |
| offset | $\delta$ | $2 \cdot 10^{-3}$ | m |
| search margin | | 0.1 | m |
| largest pressure increment | | $0.5\, D/R^3$ | Pa |

**Discretization**

The mesh is a structured mesh of quadrilaterals with $N$ cells around the ring and $L$ layers
through the wall, so that cells are square for $N t = 2 \pi R L$. The displacement is discretized
with the box method (@ref BoxDiscretization), the Jacobian is computed by numeric
differentiation, and each Newton step solves with UMFPack. The pressure is raised in increments
of at most $0.5\, D/R^3$, halved when Newton's method fails.

For each pressure the test reports the semi-minor axis of the midline, the contact force $Q$ (the
sum of the vertical contact forces on the nodes of the upper half), and the position of the
largest contact pressure along the wall. The comparison uses

- the slope $\mathrm{d}Q/\mathrm{d}p$ of a linear fit of $Q(p)$ over $6 \leq p R^3/D \leq 9$,
  where the elastica's contact force is nearly linear in $p$, and the zero of this fit, the
  extrapolated pressure of first contact,
- the onset of line contact, where the largest contact pressure leaves the minor axis.

**Results**

The table and the figures are produced in the build directory by

```bash
cd build-cmake/test/solidmechanics/contact/ring
make test_contact_ring
python3 run_benchmark.py ../../../../doc/doxygen/images
```

which runs three meshes in parallel (about 11 minutes; numpy, scipy and matplotlib) and computes
the elastica for $h/R = 0.02$.

| Mesh $N \times L$ | $\mathrm{d}Q/\mathrm{d}p$ | deviation | zero of the fit $p R^3/D$ | deviation | line contact starts at $p R^3/D$ | smallest semi-minor axis $/R$ |
| :--- | ---: | ---: | ---: | ---: | :--- | ---: |
| elastica | 1.5664 | | 5.1373 | | 10.07 | 0.02 |
| 640 × 4 | 1.5719 | +0.35 % | 5.5249 | +7.54 % | 11.5 to 12 | 0.02072 |
| 1280 × 8 | 1.5739 | +0.48 % | 5.2049 | +1.32 % | 10.5 to 11 | 0.02072 |
| 2560 × 16 | 1.5744 | +0.51 % | 5.1249 | −0.24 % | 10.5 to 11 | 0.02072 |

![Left: contact force Q between the two halves of the ring against the pressure p, for the elastica (solid line, the diamond marks the onset of line contact) and three meshes (circles dashed, squares dotted, triangles dash-dotted). Right: the difference of the computed contact force from the elastica.](ring_contactforce.png)

![The deformed ring on the finest mesh at three pressures: point contact at the lowest and line contact at the highest pressure; dots mark the nodes in contact.](ring_shapes.png)

- The extrapolated pressure of first contact converges to the elastica: it deviates by
  +7.54 %, +1.32 % and −0.24 % on the three meshes.
- The slope $\mathrm{d}Q/\mathrm{d}p$ deviates from the elastica by 0.48 % and 0.51 % on the two
  finer meshes, so the difference is not a discretization error. The elastica neglects the
  extension and the shear of the wall and its compression in the contact.
- The midline stays at least $0.02072\,R$ from the axis on all meshes, above $h = 0.02\,R$: the
  walls do not cross.
- The largest contact pressure leaves the minor axis between $p R^3/D = 10.5$ and 11 on the two
  finer meshes, after the elastica's onset of line contact at 10.07. The criteria differ: the
  elastica's line contact starts where the curvature of the midline at the contact point vanishes.

The CTest `test_contact_ring` (`ringtest.py`) runs the mesh $1280 \times 8$ for
$p R^3/D = 6, 7, \dots, 11$ and checks that the slope is within 1 % and the zero of the fit within
3 % of the elastica, that the largest contact pressure leaves the minor axis between
$p R^3/D = 10$ and 11, and that the semi-minor axis of the midline stays above $h$. Without contact
($k = 0$) the walls cross and the test fails.
