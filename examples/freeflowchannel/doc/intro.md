This example is based on dumux/test/freeflow/navierstokes/channel/2d.

# Freeflow through a channel

## Problem set-up
This example contains a stationary free flow of a fluid through two parallel solid plates in two dimensions from left to right. The figure below shows the simulation set-up. The fluid flows into the system at the left with a constant velocity of $` v = 1~\frac{\text{m}}{\text{s}} `$. The inflow velocity profile is a block profile. Due to the no-slip, no-flow boundary conditions at the top and bottom plate, the velocity profile gradually assumes a parabolic shape along the channel. At the outlet, the pressure is fixed and a zero velocity gradient in x-direction is assumed. The physical domain, which is modeled is the rectangular domain $`x\in[0,10],~y\in[0,1]`$.

![](./img/setup.png)

## Model description
The Stokes model without gravitation and without sources or sinks for a stationary, incompressible, laminar, single phase, one-component, isothermal ($`T=10^\circ C`$) flow is considered assuming a Newtonian fluid of constant density $` \varrho = 1~\frac{\text{kg}}{\text{m}^3} `$ and constant kinematic viscosity $` \nu = 1~\frac{\text{m}^2}{\text{s}} `$. The momentum balance
```math
- \nabla\cdot\left(\mu\left(\nabla\boldsymbol{u}+\nabla\boldsymbol{u}^{\text{T}}\right)\right)+ \nabla p = 0
```
with density  $`\varrho`$, velocity $`\boldsymbol{u}`$, dynamic viscosity  $`\mu=\varrho\nu`$ and pressure $`p`$ and the mass balance
```math
\nabla \cdot \left(\varrho\boldsymbol{u}\right) =0
```
are discretized using a staggered-grid finite-volume scheme as spatial discretization with pressures and velocity components as primary variables. For details on the discretization scheme, have a look at the dumux [handbook](https://dumux.org/handbook).

In the following, we take a close look at the files containing the set-up: At first, boundary conditions are set in `problem.hh` for the Navier-Stokes model. Afterwards, we show the different steps for solving the model in the source file `main.cc`. At the end, we show some simulation results.
