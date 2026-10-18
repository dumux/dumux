# Flow in a Rough Channel {#benchmark-rough-channel}

## Normal flow with bed friction after Gauckler, Manning and Strickler

**Problem Description**

A straight channel of length $L$ and width $W$ has a bed that descends along the channel with the slope
$S$. Water enters the channel with the discharge $q$ per unit width. The bed shear stress follows the law
of Gauckler, Manning and Strickler (@ref Dumux::FrictionLawManning),

$$
\frac{\tau_b}{\rho} = \frac{g n^2}{h^{1/3}} \, u |u|,
$$

where $\rho$ is the density of water, $g$ the gravitational acceleration, $n$ Manning's coefficient,
$h$ the water depth and $u$ the depth-averaged velocity along the channel. The side walls are impermeable
and free of friction.

**Exact Solution**

In normal flow, depth and velocity are uniform and the bed friction balances gravity, $g h S = \tau_b/\rho$.
With $u = q/h$, the normal depth and the normal velocity are

$$
h_n = \left( \frac{n q}{\sqrt{S}} \right)^{3/5}, \qquad u_n = \frac{q}{h_n}.
$$

For the parameters below, $h_n = 0.8685$ m and $u_n = 1.1514$ m/s.

**Parameters**

| Parameter                       | Symbol | Value      | Unit              |
|---------------------------------|--------|------------|-------------------|
| Channel length                  | $L$    | $500$      | m                 |
| Channel width                   | $W$    | $5$        | m                 |
| Bed slope                       | $S$    | $10^{-3}$  | -                 |
| Manning's coefficient           | $n$    | $0.025$    | s/m<sup>1/3</sup> |
| Discharge per unit width        | $q$    | $1$        | m<sup>2</sup>/s   |
| Gravitational acceleration      | $g$    | $9.81$     | m/s<sup>2</sup>   |

**Setup**

The @ref ShallowWaterModels "shallow water model" is discretized with cell-centered finite volumes
(@ref CCTpfaDiscretization) and the implicit Euler method. The flux over a face is the flux of the exact
solution of the Riemann problem with hydrostatic reconstruction of the bed @cite Audusse2004, and the bed
friction enters the momentum balance as a source term. The bed elevation is constant in each cell.
The discharge $q$ is prescribed at the inflow and the normal depth $h_n$ at the outflow. Since the outside
state at the boundary face has the bed of the inside cell, the prescribed depth is shifted by the bed
difference between face and cell center, such that the free surface elevation at the face is matched.

The grid of `params.input` has 500 cells along and 5 cells across the channel. The simulation starts from
the normal flow and runs up to $t = 900$ s with time steps of up to 60 s. On this grid, the $L^1$ error of
the depth at $t = 900$ s differs from that of the steady state by $3 \cdot 10^{-4}$ of its value. The test
`test_shallowwater_roughchannel` compares the solution with a reference solution. The example
[shallowwaterfriction](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/shallowwaterfriction/README.md)
explains the implementation of this setup, there starting from a depth of 1 m and water at rest.

**Results**

The figure and the tables are produced by running
```bash
python3 plot_results.py ./test_shallowwater_roughchannel
```
in the build directory `test/freeflow/shallowwater/roughchannel` (requires NumPy and Matplotlib).
The grid is refined along the channel.

![Deviation of the water depth from the normal depth (left) and L1 errors under grid refinement (right)](shallowwater_roughchannel.png){html: width=90%}

The table lists the discrete $L^1$ errors at $t = 900$ s, $\| e \|_1 = \sum_i |V_i| \, |e_i| / \sum_i |V_i|$,
where $e_i$ is the difference to the normal flow at the center of the cell $V_i$ with area $|V_i|$.

| cells along the channel | $L^1$ error $h$ (m) | rate | $L^1$ error $u$ (m/s) | rate |
|--:|--:|--:|--:|--:|
| 250 | 3.731e-04 | - | 1.543e-03 | - |
| 500 | 1.873e-04 | 0.99 | 7.700e-04 | 1.00 |
| 1000 | 9.382e-05 | 1.00 | 3.846e-04 | 1.00 |
| 2000 | 4.696e-05 | 1.00 | 1.922e-04 | 1.00 |

Depth and velocity converge with first order. Except in the cell at the inflow, the depth is below $h_n$,
by at most $5.55 \cdot 10^{-4}$ m on 250 cells and $6.96 \cdot 10^{-5}$ m on 2000 cells, in the second cell
from the inflow. Toward the outflow it approaches the prescribed normal depth; in the cell next to the
outflow it deviates by $3.6 \cdot 10^{-6}$ m on 250 cells and $5.7 \cdot 10^{-8}$ m on 2000 cells.

The test `test_shallowwater_roughchannel_convergence` runs the grids with 250, 500 and 1000 cells and checks
that the mean convergence rates of the $L^1$ errors of depth and velocity are at least 0.95.
