# Free flow through a channel

__You learn how to__

* solve a free-flow channel problem
* set outflow boundary conditions in the free-flow context

__Results__. In this example we will obtain the following stationary velocity profile:

![](./img/velocity.png)

__Table of contents__. This description is structured as follows:

[[_TOC_]]

## Mathematical model
In this example, the Stokes model for stationary and incompressible single phase flow is considered.
Thus, the momentum balance equations

```math
- \nabla\cdot\left(\mu\left(\nabla\boldsymbol{u}+\nabla\boldsymbol{u}^{\text{T}}\right)\right)+ \nabla p = 0
```

and the mass balance

```math
\nabla \cdot \left(\boldsymbol{u}\right) =0
```

are solved, where $`\varrho`$ and $`\mu`$ are the density and viscosity of the fluid,
$`\boldsymbol{u}`$ is the fluid velocity and $`p`$ is the pressure. Here, we use constant fluid
properties with $`\varrho = 1~\frac{\text{kg}}{\text{m}^3}`$ and $`\mu = 1~\text{Pa}\text{s}`$.
Furthermore, isothermal conditions with a homogeneous temperature distribution of $`T=10^\circ C`$ are assumed.

All equations are discretized with the staggered-grid finite-volume scheme as spatial discretization
with pressures and velocity components as primary variables. For details on the discretization scheme,
have a look at the Dumux [handbook](https://dumux.org/handbook).

## Problem set-up
This example considers stationary flow of a fluid between two parallel solid plates in two dimensions.
Flow is enforced from left to right by prescribing an inflow velocity of $` v = 1~\frac{\text{m}}{\text{s}} `$
on the left boundary, while a fixed pressure of $`p = 1.1 \text{bar}`$ and a zero velocity gradient
in x-direction are prescribed on the right boundary. On the top and bottom boundaries, no-slip
conditions are applied, which cause a parabolic velocity profile to develop along the channel.
Take a look at Figure 1 for an illustration of the domain and the boundary conditions.

<figure>
    <center>
        <img src="img/setup.png" alt="Free-flow setup" width="80%"/>
        <figcaption> <b> Fig.1 </b> - Setup for the free flow problem.</figcaption>
    </center>
</figure>

# Implementation

## Folder layout and files

```
└── freeflowchannel/
    ├── CMakeLists.txt          -> build system file
    ├── main.cc                 -> main program flow
    ├── params.input            -> runtime parameters
    ├── properties.hh           -> compile time configuration
    └── problem.hh              -> boundary & initial conditions
```
