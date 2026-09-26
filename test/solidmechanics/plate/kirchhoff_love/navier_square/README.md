# Simply Supported Square Kirchhoff-Love Plate {#benchmark-navier-square}

**Problem Description**

The square $0 < x, y < L$ is simply supported on all four edges and carries a uniform
load $F$ against the deflection,

$$-D\,\nabla^4 w = F \quad \text{in } \Omega, \qquad w = 0,\quad M_{nn} = 0 \quad \text{on } \partial\Omega .$$

Navier's double sine series is the exact solution,

$$w = -\frac{16FL^4}{\pi^6 D}\sum_{m,n\ \mathrm{odd}}\frac{\sin(m\pi x/L)\sin(n\pi y/L)}{mn(m^2+n^2)^2},$$

with the centre deflection $-0.00406235\,FL^4/D$, the centre moment
$M_{xx} = -0.0479\,FL^2$ at $\nu = 0.3$ and the twisting moment
$M_{xy} = \pm 0.0325\,FL^2$ at the corners. Its jump between the two edges of a corner is the
Kirchhoff corner force $2|M_{xy}| = 0.0650\,FL^2$, which holds the corners down.

**The simply supported edge in the mixed formulation**

The divergence-free tensor $\mathbf{T} = \mathbf{M}(\boldsymbol{\theta}) - \varphi\mathbf{I} - \psi\mathbf{J}$
has the traction components $\mathbf{n}\cdot\mathbf{T}\mathbf{n} = M_{nn} - \varphi$ and
$\mathbf{s}\cdot\mathbf{T}\mathbf{n} = M_{ns} - \psi$ with $\mathbf{s} = \mathbf{J}\mathbf{n}$, so the
normal traction $\mathbf{n}\cdot\mathbf{T}\mathbf{n} = -\varphi$ imposes $M_{nn} = 0$. The test
runs two choices for the tangential direction (`Problem.Mapping`):

- `tangential`: the tangential rotation $\boldsymbol{\theta}\cdot\mathbf{s} = \partial_s w = 0$
  is prescribed, which on the edges of the square is one Cartesian component, and its
  traction is the reaction. The potential $\varphi = 0$ on the edges fixes the gauge as on a
  clamped edge, and with every edge supported the constant of $\psi$ is pinned at one
  interior node.
- `traction`: the free-edge traction $\mathbf{T}\mathbf{n} = -\varphi\mathbf{n}$, which ties
  $\psi = M_{ns}$ on the edge; $\varphi$ is pinned at one support node.

At the corner $(0,0)$ the edge $y = 0$ has $M_{ns} = M_{xy}$ and the edge $x = 0$ has
$M_{ns} = -M_{xy}$, so with `traction` a single-valued $\psi$ can satisfy both edges only if
$M_{xy} = 0$ there, which contradicts the corner force of the plate.

**Results**

$L = 1$, $E = 10^6$, $t = 0.05$, $\nu = 0.3$, $F = 1$, structured triangulations of
$n\times n$ squares. Relative deviations from the series, the corner value the worst of
the four corners:

| $n$ | `tangential` $w(L/2,L/2)$ | `tangential` $M_{xy}$ corners | `traction` $w(L/2,L/2)$ | `traction` $M_{xy}$ corners |
|---|---|---|---|---|
| 8 | $5.5\times10^{-2}$ | $1.4\times10^{-1}$ | $1.1\times10^{-1}$ | $1.6\times10^{-1}$ |
| 16 | $1.4\times10^{-2}$ | $4.9\times10^{-2}$ | $3.5\times10^{-2}$ | $4.9\times10^{-2}$ |
| 32 | $3.5\times10^{-3}$ | $1.7\times10^{-2}$ | $1.0\times10^{-2}$ | $7.3\times10^{-2}$ |
| 64 | $8.9\times10^{-4}$ | $5.3\times10^{-3}$ | $2.8\times10^{-3}$ | $8.2\times10^{-2}$ |
| 128 | $2.2\times10^{-4}$ | $1.7\times10^{-3}$ | $7.8\times10^{-4}$ | $8.4\times10^{-2}$ |

With the tangential rotation prescribed the centre deflection converges at order 2.00 and
the corner twisting moment at order 1.5 to 1.7; the moment at a corner is read from the one
or two elements that touch it. With the free-edge traction the deflection still converges,
at order 1.7 to 1.9, but the twisting moment at the corners settles about 8 % away from
the series, so the corner force is wrong.

`run_navier_square.py` runs both mappings on the five meshes and fails if the
`tangential` mapping converges at less than order 1.8 in the deflection or 1.3 in the
corner moment.

Build and run from the build directory:

```sh
cmake --build . --target test_kirchhoff_love_navier_square
ctest -V -R '^test_kirchhoff_love_navier_square$'
```
