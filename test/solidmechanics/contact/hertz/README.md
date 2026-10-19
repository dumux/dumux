# Hertz Contact of Two Cylinders {#benchmark-contact-hertz}

## Contact pressure of two elastic disks pressed together, with penalty contact within one mesh

**Problem description**

Two elastic cylinders of radius $R$ and equal material touch along a line and are pressed
together. In plane strain, the contact half-width and the pressure distribution follow from
Hertz's theory @cite Hertz1882 @cite Johnson1985. The benchmark compares the contact pressure
computed with `Contact::PenaltyContact` and `Contact::ContactAssembler` with this solution. Both
bodies are parts of one mesh and their curved surfaces are one contact surface, so the test runs
the same code as self-contact.

Each body is half a disk. The upper body is the lower half of the disk with centre $(0, R)$, the
lower body the upper half of the disk with centre $(0, -R)$; they touch at the origin. The flat
edges $y = \pm R$ lie on rigid frictionless platens that move towards each other by the approach
$\Delta$ each. By symmetry, this is the problem of two full disks compressed between two pairs of
frictionless plates. The material is neo-Hookean with the first Piola-Kirchhoff stress

$$\mathbf{P} = \mu \left( \mathbf{F} - \mathbf{F}^{-T} \right) + \lambda \ln J \, \mathbf{F}^{-T},$$

where $\mathbf{F}$ is the deformation gradient, $J = \det \mathbf{F}$, and the Lamé parameters
$\mu$ and $\lambda$ follow from the Young's modulus $E$ and the Poisson ratio $\nu$. The momentum
balance $\nabla \cdot \mathbf{P} = \mathbf{0}$ holds in the reference configuration.

![Two half disks of radius R between frictionless platens that move towards each other by the approach Delta. The dash-dotted line is the axis x = 0, on which the horizontal displacement is zero. The inset shows the contact zone of width 2a with the elliptic pressure distribution of peak value p0.](hertz_setup.svg){html: width=80%}

| Boundary | Condition |
|----------|-----------|
| flat edge $y = R$ of the upper body | $u_y = -\Delta$, no tangential traction |
| flat edge $y = -R$ of the lower body | $u_y = \Delta$, no tangential traction |
| nodes on the axis $x = 0$ | $u_x = 0$ |
| curved surfaces | traction free, except for the contact pressure |

**Penalty contact**

Every node of the contact surface is tested against the opposite facets. The gap $g$ is the signed
distance from the plane of a facet, negative behind it. The contact pressure is a continuously
differentiable function of the gap,

$$p_\mathrm{c}(g) = \begin{cases}
0 & g > \delta, \\
k (\delta - g)^2 / (2\delta) & 0 \leq g \leq \delta, \\
k (\delta/2 - g) & g < 0,
\end{cases}$$

with the penalty stiffness $k$ and the offset $\delta$, at which the contact engages. The contact
energy of a node is $\frac{1}{2} A_s \psi(g)$, with $\psi' = -p_\mathrm{c}$ and the length $A_s$ of the dual
surface patch of the node in the reference configuration. The penalty pressure is therefore a
force per reference length. Contact forces and their Jacobian are added to the system of the
hyperelastic model by `Contact::ContactAssembler`; the class documentation of
`Contact::PenaltyContact` gives the candidate search and the smoothing between facets.

**Reference solution**

For the force $F$ per unit length between the bodies, Hertz's theory gives the contact
half-width and the pressure

$$a = \sqrt{\frac{4 F R (1 - \nu^2)}{\pi E}}, \qquad p(x) = p_0 \sqrt{1 - \frac{x^2}{a^2}}, \qquad
p_0 = \frac{2F}{\pi a},$$

for $|x| \leq a$, and zero outside. The theory assumes small strains and half-spaces. The strain
at the contact is of the order $p_0/E = a / (2 R (1 - \nu^2))$, so terms of this order, about
$a/(2R)$, are beyond the accuracy of the reference.

**Parameters**

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| radius | $R$ | 1 | m |
| Young's modulus | $E$ | 1 | Pa |
| Poisson ratio | $\nu$ | 0.3 | - |
| penalty stiffness | $k$ | $10^5$ | Pa/m |
| offset | $\delta$ | $5 \cdot 10^{-4}$ | m |
| search margin | | 0.01 | m |
| largest approach increment | | $5 \cdot 10^{-4}$ | m |

**Discretization**

The meshes are unstructured triangle meshes (Gmsh, `hertz.geo`) with the mesh size $h_c$ within
$0.15 R$ of the origin, growing to $R/10$ at $0.8 R$. Gmsh places the surface nodes of the two
arcs at the same horizontal positions, so the nodes of the two bodies face each other. The
displacement is discretized with the box
method (@ref BoxDiscretization), the Jacobian is computed by numeric differentiation, and each
Newton step solves with UMFPack. The contact surface consists of the boundary facets with
$|y| < R/2$. The approach is raised in increments of at most $5 \cdot 10^{-4}\,R$, halved when
Newton's method fails.

For each approach, the test computes

- the force $F$, the sum of the vertical contact forces on the nodes of the upper body, and from it
  $a$ and $p_0$,
- the nodal pressure $p_s = p_\mathrm{c}(g_s) A_s / A_s^\mathrm{cur}$, the penalty pressure converted to a
  force per current length with the current length $A_s^\mathrm{cur}$ of the dual surface patch,
- the peak error $p_\mathrm{max}/p_0 - 1$ with the largest nodal pressure $p_\mathrm{max}$,
- the relative $L^2$ error over the surface nodes with $|x_s| < 3a$,
  $$e = \left( \frac{\sum_s A_s^\mathrm{cur} \left( p_s - p(x_s) \right)^2}{\sum_s A_s^\mathrm{cur} p(x_s)^2} \right)^{1/2},$$
  where $x_s$ is the current horizontal position of node $s$.

**Results**

The tables and the figures are produced in the build directory by

```bash
cd build-cmake/test/solidmechanics/contact/hertz
make test_contact_hertz
python3 run_benchmark.py $(which gmsh)
```

which meshes the geometry with the mesh sizes of the tables and runs them in parallel (about 16
minutes; numpy and matplotlib); with the directory `doc/doxygen/images` of the DuMux sources as
second argument, it also writes the figures. The approaches $\Delta/R = 0.00125$, 0.0025, 0.005 and
0.01 give contact half-widths $a/R$ from 0.026 to 0.077. The force $F/(ER)$ agrees between the
meshes to within 0.04 %.

Peak pressure, $p_\mathrm{max}/p_0 - 1$:

| $a/R$ | $p_0/E$ | $h_c = 0.005\,R$ | $h_c = 0.0025\,R$ | $h_c = 0.00125\,R$ |
| ---: | ---: | ---: | ---: | ---: |
| 0.0257 | 0.0141 | −0.16 % | +0.17 % | +0.17 % |
| 0.0361 | 0.0199 | +0.29 % | +0.54 % | +0.56 % |
| 0.0523 | 0.0287 | +0.69 % | +0.94 % | +1.01 % |
| 0.0771 | 0.0423 | +1.64 % | +1.52 % | +1.60 % |

Relative $L^2$ error of the pressure profile:

| $a/R$ | $p_0/E$ | $h_c = 0.005\,R$ | $h_c = 0.0025\,R$ | $h_c = 0.00125\,R$ |
| ---: | ---: | ---: | ---: | ---: |
| 0.0257 | 0.0141 | 0.75 % | 1.10 % | 0.79 % |
| 0.0361 | 0.0199 | 0.80 % | 0.86 % | 0.76 % |
| 0.0523 | 0.0287 | 1.86 % | 1.88 % | 1.53 % |
| 0.0771 | 0.0423 | 3.03 % | 2.66 % | 2.46 % |

![Contact pressure over the current horizontal position on the finest mesh for a/R = 0.036 and 0.077 (circles and squares) and Hertz's solution (solid line), both scaled with a and p0 of the computed force.](hertz_profiles.png){html: width=65%}

![Peak error of the nodal pressure for the mesh sizes of the tables (circles solid, squares dashed, triangles dotted) and of the penalty pressure, a force per reference length, on the finest mesh (crosses dash-dotted), and the relative L2 error of the pressure profile.](hertz_errors.png)

- The computed pressure follows Hertz's ellipse (pressure profiles above). The deviations do not
  decrease under mesh refinement: the peak error changes by at most 0.08 percentage points between
  the meshes with $h_c = 0.0025\,R$ and $0.00125\,R$, and the $L^2$ error stays between 0.8 % and
  2.5 % on the finest mesh.
- The deviations grow with the load. They are of the size of the terms Hertz's small-strain solution
  neglects: the peak of the penalty pressure, which is the same pressure per reference instead of
  current length, lies 0.43 % to 0.54 % below $p_0$, and differs from the peak of the nodal pressure
  by 0.71 %, 1.00 %, 1.44 % and 2.09 %, which is $p_0/(2E)$ to within 0.03 percentage points.
  Hertz's linear theory does not distinguish the two, and $p_0$ lies between them for every load on
  the meshes with $h_c = 0.0025\,R$ and $0.00125\,R$.

The CTest `test_contact_hertz` (`hertztest.py`) meshes with $h_c = 0.005\,R$, runs the approaches
0.0025 and 0.005, and checks that the peak pressure is within 2 % of $p_0$ and that the $L^2$
error is below $a/R$.
