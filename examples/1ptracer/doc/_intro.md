# Single-phase flow and tracer transport

In this example, single-phase flow and tracer transport through a domain with a
heterogeneous permeability distribution is considered. A velocity distribution
is obtained from the solution of a stationary single-phase problem, and subsequently,
this velocity field is used for the simulation of the transport of tracer through the
domain.

__The main points illustrated in this example are__
* setting up and solving a stationary single-phase flow problem
* setting up and solving a tracer transport problem
* solving two problems sequentially and realizing the data transfer
* using a simple method to generate a random permeability field

__Table of contents__. This description is structured as follows:

[[_TOC_]]

## Problem set-up
A domain with an extent of $`10 \, \mathrm{m} \times 10 \, \mathrm{m}`$ is considered,
in which a heterogeneous permeability distribution is generated randomly. The problem
set-up is shown in the figure below. In the stationary single-phase simulation,
a pressure difference between the bottom and the top boundaries is prescribed (see left figure),
resulting in a non-uniform velocity distribution due to the heterogeneous medium.
Neumann no-flow conditions are used on the lateral sides. On the basis of the resulting
velocity field, the transport of an initial tracer concentration distribution through
the domain is simulated (see right figure). Initially, non-zero tracer concentrations
are prescribed on a small strip close to the bottom boundary.

<figure>
    <center>
        <img src="img/1p_setup.png" alt="Single-phase setup" width="47%"/> <img src="img/xtracer.gif" alt="Tracer result" width="47%"/>
        <figcaption> <b> Fig.1 </b> - Setup for the single-phase problem (left) and tracer mass fraction over time as computed with the tracer model (right).</figcaption>
    </center>
</figure>


## Model description

As mentioned above, two models are solved sequentially in this example. A single-phase
model (_1p model_) is used to solve for the stationary velocity distribution of a fluid phase
in the domain. The tracer transport is solved with the _tracer model_, which solves an advection-diffusion
equation for a tracer component, which is assumed not to affect the density and viscosity
of the fluid phase.

### 1p Model
The single phase model uses Darcy's law as the equation for the momentum conservation:

```math
\textbf v = - \frac{\textbf K}{\mu} \left(\textbf{grad}\, p - \varrho {\textbf g} \right),
```

with the darcy velocity $`\textbf v`$, the permeability $`\textbf K`$, the dynamic viscosity $`\mu`$, the pressure $`p`$, the density $`\varrho`$ and the gravitational acceleration $`\textbf g`$.

Darcy's law is inserted into the mass balance equation:

```math
\phi \frac{\partial \varrho}{\partial t} + \text{div} \left( \varrho \textbf v \right) = 0,
```

where $`\phi`$ is the porosity. The primary variable used in this model is the pressure $`p`$.

### Tracer Model
The tracer model solves the mass conservation equation of a tracer component $`\kappa`$,
in which both advective and diffusive transport mechanisms are considered:

```math
\phi \frac{ \partial \varrho X^\kappa}{\partial t} - \text{div} \left\lbrace \varrho X^\kappa {\textbf v} + \varrho D^\kappa_\text{pm} \textbf{grad} X^\kappa \right\rbrace = 0.
```

Here, $`\textbf v`$ is a velocity field, which in this example is computed using the _1p model_ (see above). Moreover, $`X^\kappa`$ is the tracer mass fraction and $` D^\kappa_\text{pm} `$ is the
effective diffusivity. In this example, the effective diffusivity is a function of the diffusion
coefficient of the tracer component $`D^\kappa`$ and the porosity and tortuosity $`\tau`$ of the porous
medium (see [dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh)):

```math
D^\kappa_\text{pm}= \phi \tau D^\kappa.
```

The primary variable used in this model is the tracer mass fraction $`X^\kappa`$.

### Discretization

In this example, all equations are discretized using cell-centered finite volumes with two-point flux
approximation as spatial discretization scheme. For details on the discretization schemes available in
DuMuX, have a look at the [handbook](https://dumux.org/handbook). We use the implicit Euler method as
time discretization scheme for the tracer component balance equation solved in the _tracer model_.

# Implementation

In the following, we take a closer look at the source files for this example:

```
└── 1ptracer/
    ├── CMakeLists.txt          -> build system file
    ├── main.cc                 -> main program flow
    ├── params.input            -> runtime parameters
    ├── properties._1p.hh       -> compile time settings for the single-phase flow simulation
    ├── problem_1p.hh           -> boundary & initial conditions for the single-phase flow simulation
    ├── spatialparams_1p.hh     -> parameter distributions for the single-phase flow simulation
    ├── properties_tracer.hh    -> compile time settings for the tracer transport simulation
    ├── problem_tracer.hh       -> boundary & initial conditions for the tracer transport simulation
    └── spatialparams_tracer.hh -> parameter distributions for the tracer transport simulation
```

In order to define a simulation setup in DuMuX, you need to implement compile-time settings,
where you specify the classes and compile-time options that DuMuX should use for the simulation.
Moreover, a `Problem` class needs to be implemented, in which the initial and boundary conditions
are specified. Finally, spatially-distributed values for the parameters required by the used model
are implemented in a `SpatialParams` class.

In the documentations behind the links provided in the following, you will find how the above-
mentioned settings and classes are realized in this example. Finally, it is discussed how the two
simulations are solved sequentially in the main file (`main.cc`).
