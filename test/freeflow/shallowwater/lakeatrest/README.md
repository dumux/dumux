# Lake at Rest {#benchmark-lake-at-rest}

## Well-balancing of the bed slope source term

**Problem Description**

Water at rest over a bed of variable elevation $z(x, y)$ is a steady state of the shallow water equations:
for a horizontal free surface $h + z = H$ and zero velocity, the gradient of the hydrostatic pressure term
$\frac{1}{2} g h^2$ in the momentum flux is balanced by the bed slope source term $-g h \nabla z$.
Here $h$ is the water depth, $H$ the elevation of the free surface and $g$ the gravitational acceleration.
A scheme that preserves this state up to round-off satisfies the C-property @cite Bermudez1994 and is
called well-balanced. A scheme that is not well-balanced generates velocities of the order of the
bed variation.

The bed on the domain $[0, 2\,\mathrm{m}] \times [0, 1\,\mathrm{m}]$ consists of a smooth hump, the bed of
the still water test of LeVeque @cite LeVeque1998, and a sill with discontinuous flanks,

$$
z(x, y) = a \exp\left( -5\,\mathrm{m}^{-2} (x - 0.9\,\mathrm{m})^2 - 50\,\mathrm{m}^{-2} (y - 0.5\,\mathrm{m})^2 \right)
        + \begin{cases} b, & 1.4\,\mathrm{m} < x < 1.6\,\mathrm{m}, \\ 0, & \text{otherwise,} \end{cases}
$$

with the hump height $a$ and the sill height $b$. All boundaries are walls and there is no friction.

![Setup of the lake at rest: side view of the bed along y = 0.5 m with both water levels, and top view with the contours at 50 % (solid) and 10 % (dashed) of the hump height; the vertical scales differ](shallowwater_lakeatrest_setup.svg){html: width=70%}

**Exact Solution**

The exact solution is the initial state, $h = \max(H - z, 0)$ and $u = v = 0$, where $u$ and $v$ are the
velocity components in $x$- and $y$-direction. For $H = 1$ m the bed is submerged everywhere. For
$H = 0.5$ m the top of the hump and the sill fall dry.

**Parameters**

| Parameter                         | Symbol           | Value                | Unit             |
|-----------------------------------|------------------|----------------------|------------------|
| Hump height                       | $a$              | $0.8$                | m                |
| Sill height                       | $b$              | $0.6$                | m                |
| Free surface elevation            | $H$              | $1$ (submerged), $0.5$ (emerged) | m    |
| Water film on the dry bed         | $h_\mathrm{min}$ | $0$ (submerged), $10^{-5}$ (emerged) | m |
| Gravitational acceleration        | $g$              | $9.81$               | m/s<sup>2</sup>  |
| Final time                        | $t_\mathrm{end}$ | $100$                | s                |
| Time step size                    | $\Delta t$       | $10$                 | s                |

**Setup**

The @ref ShallowWaterModels "shallow water model" is discretized with cell-centered finite volumes
(@ref CCTpfaDiscretization) on a grid of $100 \times 50$ cells and the implicit Euler method. The bed
elevation is constant in each cell and is evaluated at the cell center. The flux over a face is computed
with the hydrostatic reconstruction of Audusse et al. @cite Audusse2004. The depths on both sides are reduced
to the free surface above the higher of the two beds, $h^* = \max(0, h + z - \max(z_\mathrm{in}, z_\mathrm{out}))$,
the flux of the exact solution of the Riemann problem between the reconstructed states is evaluated,
and the momentum flux is corrected by $\frac{1}{2} g (h^{*2} - h^2)$. At rest, the reconstructed states on
both sides of a face are equal, and the corrected momentum flux out of a cell is $\frac{1}{2} g h^2$ of that
cell on each of its faces, which sums to zero over the cell, also across the flanks of the sill.

A completely dry cell has no momentum storage, which makes the Jacobian of the momentum balance singular.
Where the bed emerges, the depth is therefore set to a film of thickness $h_\mathrm{min}$
(`Problem.MinWaterDepth`). The film on the sloping bed is not at rest itself, which limits the accuracy that
can be reached on a partially dry bed.

The test `test_shallowwater_lakeatrest` (and `_parallel` on two processes) runs the submerged case and fails
if the maximum deviation of the free surface elevation or of the specific discharge $h \sqrt{u^2 + v^2}$
from the exact solution exceeds $10^{-10}$ at $t_\mathrm{end}$. The test
`test_shallowwater_lakeatrest_emerged` runs the emerged case with the tolerance $10^{-5}$.

**Results**

The figure and the table are produced by running
```bash
python3 plot_results.py ./test_shallowwater_lakeatrest
```
in the build directory `test/freeflow/shallowwater/lakeatrest` (requires NumPy and Matplotlib).

![Bed and free surface along y = 0.49 m, and deviation of the water depth from the exact solution relative to the water level, at t = 100 s](shallowwater_lakeatrest.png){html: width=90%}

The table lists the maximum deviations at $t_\mathrm{end}$ relative to the scales of the problem: the
deviation of the free surface elevation relative to $H$, and the deviation of the specific discharge
relative to $H \sqrt{g H}$, the discharge of water of depth $H$ moving with the celerity $\sqrt{g H}$.

| case | $H$ (m) | $h_\mathrm{min}$ (m) | max. deviation free surface / $H$ | max. deviation specific discharge / $H \sqrt{g H}$ |
|:--|--:|--:|--:|--:|
| submerged | 1.0 | 0 | 9.4e-15 | 1.0e-14 |
| emerged | 0.5 | 1e-05 | 4.4e-07 | 3.9e-07 |

On the submerged bed the lake at rest is preserved to round-off. On the partially dry bed the deviation
of the free surface is $2.2 \cdot 10^{-7}$ m, 2.2 % of the film thickness.
