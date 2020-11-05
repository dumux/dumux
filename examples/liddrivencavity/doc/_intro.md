# Shear-driven cavity flow

We use the Navier-Stokes equations to simulate laminar incompressible flow in a
cavity whose lid moves with a constant velocity u = 1 m/s.
We will verify the numerical model by comparing the simulation results with reference data published in [Ghia et al. (1982)](https://doi.org/10.1016/0021-9991(82)90058-4) and [Jurjević (1999)](https://doi.org/10.1002/(SICI)1097-0363(19991015)31:3<601::AID-FLD892>3.0.CO;2-Z).

__Results__. After simulating a few time steps, we will obtain the following velocity field for Reynolds number Re = 1 and Re = 1000:
<figure>
    <center>
        <img src="img/result.svg" alt="Numerical results" width="50%"/>
        <figcaption> <b> Fig.1 </b> - Steady velocity field for Stokes (left) and Navier-Stokes flow for the lid-driven cavity problem.</figcaption>
    </center>
</figure>
Our numerical results agree well with the reference data:
<figure>
    <center>
        <img src="img/lidverification.png" alt="Lid-driven cavity verification" width="90%"/>
    </center>
    <figcaption> <b> Fig.2 </b> - Horizontal and vertical velocity profiles at x = 0.5 m and y = 0.5 m for Re = 1 (left) and Re = 1000 (right). exp: experimental data; num: numerical data.</figcaption>
</figure>

__In this example, you will__

* solve a single-phase Navier-Stokes flow problem
* see the differences between Stokes flow (Re = 1) and Navier-Stokes flow (Re = 1000)
* compare the numerical results with the reference data using the plotting library `matplotlib`

__Table of contents__. This description is structured as follows:

[[_TOC_]]

## Problem setup

Flow in a cavity with the dimensions 1 m × 1 m is considered, where the top lid is
moving at a constant speed of 1 m/s to the right.

The following figure illustrates the setup:

<figure>
    <center>
        <img src="img/setup.png" alt="Lid-driven cavity setup" width="35%"/>
        <figcaption> <b> Fig.3 </b> - Setup for the lid-driven cavity problem.</figcaption>
    </center>
</figure>

Two different flow regimes at Re = 1 ($`\nu = 1`$ $`\text{m}^2/\text{s}`$) and Re = 1000 ($`\nu = 1000`$  $`\text{m}^2/\text{s}`$) are simulated, where the Reynolds number is defined with respect to the cavity’s side length.

## Mathematical & numerical models

Mass and momentum balance are given by

```math
\nabla \cdot \bold{v} =0,
```
```math
 \frac{(\partial\rho\bold{v})}{\partial t} + \nabla \cdot (\rho\bold{v}\bold{v}^{\text{T}}) =\nabla\cdot\left[\mu\left(\nabla\bold{v}+\nabla\bold{v}^{\text{T}}\right)\right]- \nabla p,
```

where $`\bold{v}`$ and p are the velocity and pressure of the fluid (primary variables). $`\rho`$ and $`\mu=\rho\nu`$ are the mass density and dynamic viscosity (fluid properties).

All equations are discretized with the staggered-grid finite-volume scheme as spatial discretization with pressures and velocity components as primary variables. For details on the discretization scheme, we refer to the DuMu<sup>x</sup> [handbook](https://dumux.org/docs/handbook/master/dumux-handbook.pdf).

# Implementation & Postprocessing
