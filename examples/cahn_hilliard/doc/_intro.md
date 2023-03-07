# Cahn-Hilliard model

A random initial distribution of two phases separating according to the Cahn-Hilliard model.

__The main points illustrated in this example are__
* A base setup for solving a nonlinear partial differential equation

__Table of contents__. This description is structured as follows:

[[_TOC_]]

__Result__. The result will look like below.
You see the concentration variable $c$ capturing the phase distribution after 1 second.

<figure>
    <center>
        <img src="img/results_phase_distribution.png" alt="Concentration field modeling phase distribution" width="50%"/>
    </center>
</figure>

## Problem set-up

A square two-dimensional domain is initialized with a concentration `c`, randomly perturbed by
noise following a uniform distribution between 0.41 and 0.43. 
With time the concentration field evolves towards attaining mostly values near to 0 or 1 while
conserving the total concentration, modeling the separation of two immiscible fluids.

## Model description

The Cahn-Hilliard model uses a pair of second order nonlinear partial differential equations,
for a concentration field $c$ and a chemical potential $\mu$. The former is conserved and
evolves with a flux proportional to $\nabla \mu$, while the latter depends on the laplacian of
the concentration $\nabla^2 c$ and a nonlinear function of $c$, a derivative of a free energy
functional.

```math
\frac{\partial c}{\partial t} = M \nabla^2 \mu \\
\mu = E (4 c^3 - 6 c^2 +2 c) - \gamma \nabla^2 c
```

Here $M$ denotes the mobility of the concentration, while the energy scale $E$ and surface tension
$\gamma$ balance the two contributions of the free energy functional.

## Implementation of a simple nonlinear PDE

For this example the C++ code is contained in two files, `main.cc` and `model.hh`.
The `model.hh` header contains the `ModelTraits` and properties for the model TypeTag
`CahnHilliardModel`, as well as volumevariables `CahnHilliardModelVolumeVariables`
and the basic local residual `CahnHilliardModelLocalResidual`.
The residual's storage term is extended in the problem definition `CahnHilliardTestProblem`
in the file `main.cc`, which also contains more specific properties of the problem setup (for
type tag `CahnHilliardTest`) as well as the actual simulation loop usually found in a mainfile.

More details are given in [main.cc](doc/mainfile.md) and [model.hh](doc/modelheader.md).

Additionally the folder contains `params.input` containing runtime parameters and the build
system file `CMakeLists.txt` defining conditions for the accompanying test.
