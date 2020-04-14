# Shallow water flow with bottom friction
In this example, the shallow water flow model is applied to simulate
a steady subcritical flow including bottom friction (bed shear stress).

__You will learn how to__

* solve a shallow water flow problem including bottom friction
* compute and output (VTK) an analytical reference solution

__Result__. The numerical and analytical solutions for the free surface will look like this:

<figure>
    <center>
        <img src="img/swe_result.png" alt="Shallow water result" width="60%"/>
        <figcaption> <b> Fig.1 </b> - Setup and result for the shallow water problem with bottom friction.</figcaption>
    </center>
</figure>

__Table of contents__. This description is structured as follows:

[[_TOC_]]

## Problem set-up
### Model domain
The model domain is given by a rough channel with a slope of 0.001.
The domain is 500 meters long and 5 meters wide.
The bottom altitude is 10 m at the inflow and hence 9.5 m at the outflow.
Bottom friction is considered by applying
Manning's law ($`n`$ = 0.025).

### Boundary conditions
At the lateral sides a no-flow boundary condition is applied. Also no friction is
considered there and therefore a no slip boundary
condition is applied. These are the default boundary condition for the shallow
water model. At the left border a discharge boundary condition
is applied as inflow boundary condition with $`q = -1.0 m^2 s^{-1}`$.
At the right border a fixed water depth boundary condition
is applied for the outflow. Normal flow is assumed, therefore the water
depth at the right border is calculated using the equation
of Gauckler, Manning and Strickler.

### Initial conditons
The initial water depth is set to 1 m, which is slightly higher than the normal flow
water depth (0.87 m). Therefore, we expect a decreasing
water level during the simulation until the normal flow condition is reached in
the entire model domain. The inital velocity is set to zero.

## Model description
As mentioned above, this examples uses the shallow water equations (SWEs) to solve the problem.
These are a depth averaged simplification of the Navier-Stokes equations. To calculate the
bottom friction Manning's law is used. An alternative is Nikuradse's law, which is also implemented
in DuMu<sup>x</sup>.

### Shallow water model
The shallow water equations are given as:

```math
\frac{\partial \mathbf{U}}{\partial t} +
\frac{\partial \mathbf{F}}{\partial x} +
\frac{\partial \mathbf{G}}{\partial y} - \mathbf{S_b} - \mathbf{S_f} = 0
```

where $`\mathbf{U}`$, $`\mathbf{F}`$ and $`\mathbf{G}`$ defined as

```math
\mathbf{U} = \begin{bmatrix} h \\ uh \\ vh \end{bmatrix},
\mathbf{F} = \begin{bmatrix} hu \\ hu^2  + \frac{1}{2} gh^2 \\ huv \end{bmatrix},
\mathbf{G} = \begin{bmatrix} hv \\ huv \\ hv^2  + \frac{1}{2} gh^2 \end{bmatrix}
```

$`h`$ the water depth, $`u`$ the velocity in x-direction and $`v`$ the velocity in y-direction,
$`g`$ is the constant of gravity.

The source terms for the bed slope $`\mathbf{S_b}`$ and friction
$`\mathbf{S_f}`$ are given as

```math
\mathbf{S_b} = \begin{bmatrix} 0 \\ -gh \frac{\partial z}{\partial x}
               \\ -gh \frac{\partial z}{\partial y}\end{bmatrix},
\mathbf{S_f} = \begin{bmatrix} 0 \\ghS_{fx} \\ghS_{fy}\end{bmatrix}.
```

with the bedSurface $`z`$. $`S_{fx}`$ and $`S_{fy}`$ are the bed shear stess
components in x- and y-direction, which are calculated by Manning's law.

### Mannings law
The empirical Manning model specifies the bed shear stress by the following equations:

```math
S_{fx} = \frac{n^2u}{R_{hy}^{4/3}} \sqrt(u^2 + v^2),

S_{fy} = \frac{n^2v}{R_{hy}^{4/3}} \sqrt(u^2 + v^2)
```

$`n`$ is Manning's friction value and $`R_{hy}`$ is the hydraulic radius,
which is assumed to be equal to the water depth $`h`$.

### Analytical solution
Since normal flow conditions are assumed, the analytic solution is calculated using the equation
of Gauckler, Manning and Strickler:

```math
v_m = n^{-1} R_{hy}^{2/3} I_s^{1/2}
```

Where the mean velocity $`v_m`$ is given as

```math
v_m = \frac{q}{h}
```

$`I_s`$ is the bed slope and $`q`$ the unity inflow discharge.

Hence, the water depth $`h`$ can be calculated by

```math
h = \left(\frac{n q}{\sqrt{I_s}} \right)^{3/5}
```

### Discretisation
For this example, a cell-centered finite volume method (cctpfa) is applied to solve the SWEs
in combination with a fully-implicit time discretization. For cases where no sharp fronts or
traveling waves occur it is possible to apply time steps larger than CFL number = 1 to reduce
the computation time. Even if a steady state solution is considered, an implicit time stepping method
is applied.

# Implementation

## Folder layout and files

```
└── shallowwaterfriction/
    ├── CMakeLists.txt          -> build system file
    ├── main.cc                 -> main program flow
    ├── params.input            -> runtime parameters
    ├── properties.hh           -> compile time configuration
    ├── problem.hh              -> boundary & initial conditions
    └── spatialparams.hh        -> spatial parameter fields
```
