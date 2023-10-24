# Cahn-Hilliard equation

In this example, we create an application solving the Cahn-Hilliard equation which describes
a phase separation process.
This example is implemented in three files, [`model.hh`](model.hh), [`params.input`](params.input), [`main.cc`](main.cc).
The [`model.hh`](model.hh) header implements a simple self-contained Cahn-Hilliard equation model
to be used with a control-volume finite element discretization like the Box method.
The [`main.cc`](main.cc) source file contains the main program and [`params.input`](params.input) contains
runtime parameters. The executable is configured in [`CMakeLists.txt`](CMakeLists.txt) and created with CMake.
The source code will be discussed in Part 1 and Part 2 after an introduction to the equations and the problem description.

__The main points illustrated in this example are__

* setting up a model with two coupled nonlinear equations
* creating a vector-valued random initial solution in parallel
* solving a time-dependent nonlinear problem
* using a Newton solver instance
* enabling partial reassembly for the Newton method
* using adaptive time stepping

__Results__. The result will look like below.

<figure>
    <center>
        <img src="img/animation.gif" alt="Cahn Hilliard result" width="400px"/>
        <figcaption> <b> Fig.1 </b> - Evolution of an initially random concentration field according to the Cahn-Hilliard equation.</figcaption>
    </center>
</figure>

__Table of contents__. This description is structured as follows:

[TOC]

## Equation and problem description

The Cahn-Hilliard equation is a forth order partial differential equation (PDE) and reads

```math
\frac{\partial c}{\partial t} - \nabla\cdot \left( M \nabla\left[ \frac{\partial f}{\partial c} - \gamma \nabla^2 c \right]\right) = 0 \quad \text{in} \; \Omega \times (0, T]
```

with the concentration $c(x,t)$, the mobility coefficient $M$, and surface tension $\gamma$.
The domain $\Omega \subset \mathbb{R}^2$ is initialized with a concentration field
$c(x,t=0) = 0.42 + \zeta$, randomly perturbed by
noise $\zeta$ following a uniform distribution $\zeta \sim U(-0.02, 0.02)$.
Over time, the concentration field evolves towards attaining mostly values near to $0$ or $1$ while
conserving the total concentration. The model describes the separation of two immiscible fluids.

The fourth order PDE cannot be solved by a standard finite volume scheme. We therefore
split the equation in two second order PDEs

```math
\begin{aligned}
\frac{\partial c}{\partial t} - \nabla\cdot \left( M \nabla \mu \right) &= 0\\
- \gamma \nabla^2 c &= \mu - \frac{\partial f}{\partial c}.
\end{aligned}
```

were $\mu$ is called chemical potential.
Using the free energy functional $f = E c^2(1-c)^2$, where $E$ is a scaling parameter,
we obtain

```math
\begin{aligned}
\frac{\partial c}{\partial t} - \nabla\cdot \left( M \nabla \mu \right) &= 0\\
- \gamma \nabla^2 c &= \mu - E (4 c^3 - 6 c^2 +2 c).
\end{aligned}
```

The concentration field is conserved and
evolves with a flux proportional to $\nabla \mu$, while the latter depends on the Laplacian of
the concentration $\nabla^2 c$ and the free energy functional which is a nonlinear function of $c$. We therefore have a system of two coupled nonlinear PDEs. We will use homogeneous Neumann
boundary conditions everywhere on $\partial \Omega$. This means that the total concentration in $\Omega$ stays constant.

For the implementation, it is useful to write the conservation equation in the following abstract
form in vector notation,

```math
\frac{\partial \boldsymbol{S}(\boldsymbol{U})}{\partial t} + \nabla\cdot \boldsymbol{F}(\boldsymbol{U}) = \boldsymbol{Q}(\boldsymbol{U})
```

where for the Cahn-Hillard equation, we have

```math
\boldsymbol{U} := \begin{bmatrix} c \\ \mu \end{bmatrix}, \quad
\boldsymbol{S}(\boldsymbol{U}) := \begin{bmatrix} c \\ 0 \end{bmatrix}, \quad
\boldsymbol{F}(\boldsymbol{U}) := \begin{bmatrix} -M \nabla \mu \\ -\gamma \nabla c \end{bmatrix}, \quad
\boldsymbol{Q}(\boldsymbol{U}) := \begin{bmatrix} 0 \\ \mu - E (4 c^3 - 6 c^2 +2 c) \end{bmatrix}.
```

### Model parameters

* $\Omega = [0,1]\times[0,1]$ (unit square)
* End time $T = 1$, initial time step size $\Delta t = 0.0025$, and maximum time step size $\Delta t = 0.01$
* $100 \times 100$ cells (structured Cartesian grid)
* $M = 0.0001$
* $\gamma = 0.01$
* $f = E c^2(1-c)^2$, with $E = 100$
* $c^0_{h,B} = 0.42 + \zeta$, where $\zeta \sim U(-0.02,0.02)$
* $\mu^0_{h,B} = 0.0$
* homogeneous Neumann boundary conditions

### Running the example in parallel

See [Diffusion equation example](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/diffusion/README.md).

# Implementation

For this example the C++ code is contained in two files, [`model.hh`](model.hh) and [`main.cc`](main.cc). The [`model.hh`](model.hh) header contains the `ModelTraits` and properties for the
model type tag `CahnHilliardModel`, as well as the volume variables class
`CahnHilliardModelVolumeVariables` and the local residual class `CahnHilliardModelLocalResidual`.
The source term is extended in the problem definition `CahnHilliardTestProblem`
in the file `main.cc`, which also contains more specific properties of the problem setup (for
type tag `CahnHilliardTest`) as well as the actual simulation loop usually found in the main function.
