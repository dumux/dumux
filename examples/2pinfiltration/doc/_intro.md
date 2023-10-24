# Two-phase flow infiltration problem with adaptive grid

__In this example, you will learn how to__

* solve a two-phase flow in porous media problem with two immiscible phases
* set boundary conditions and a simple injection well
* specify a point source
* read the initial solution from a text file
* specify a lens with different porous material parameters
* use adaptive grid refinement around the saturation front


__Prerequisites__. You need [dune-alugrid](https://gitlab.dune-project.org/extensions/dune-alugrid) in order to compile and run this example.

__Result__. The resulting saturation distribution in this example will look like this:

![](./img/test_2p_pointsource_adaptive.png)

__Table of contents__. This description is structured as follows:

[[_TOC_]]

## Scenario and mathematical model

We model a soil contamination problem where DNAPL infiltrates a porous medium. The initial distribution of DNAPL is known and we can read it from a txt-file.
We describe the problem using a two phase model with two immiscible fluid phases (subscripts $`w`$ and $`n`$). We use multiphase Darcy's law for the momentum balance equations of the fluid phases, i.e.:

```math
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
```

Inserting this into the conservation equations for each phase $`\alpha`$ yields:

```math
\phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -\textbf{div} \left\{ \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \bf g \right)
 \right\} - q_\alpha = 0
```

To reduce the number of unknowns and close the system we need closure relations for these equations. In this example, we use the
Van Genuchten-Mualem relationships (see
[Van Genuchten (1980)](https://doi.org/10.2136/sssaj1980.03615995004400050002x)
and
[Mualem (1976)](https://doi.org/10.1029/WR012i003p00513))
for the capillary pressure $`pc = p_n - p_w`$ and the relative permeabilities $`k_r\alpha`$.
The parameters for these relationships are specified in the `spatialparams.hh` file.

With the additional constraint that $`S_w + S_n = 1`$, the number of unknowns is reduced to two.
In this example we use the wetting phase pressure $`p_w`$ and the saturation of the nonwetting phase $`S_n`$ as primary variables. It is also possible to switch that formulation to the nonwetting pressure and the wetting saturation.

The two-dimensional model domain is 6m x 4m and contains a lens with a lower permeability and porosity. We read the initial values for the DNAPL saturation and the water pressure from a file.
The lens and the initial saturation can be seen in Figures 1 and 2.

![](./img/test_2p_pointsource_lens.png)

At the left and the right boundary of the domain we set a linear pressure gradient as a Dirichlet boundary condition. On the upper and lower boundary we set Neumann boundary conditions.
DNAPL enters the model domain at the upper boundary between 1.75m ≤ x ≤ 2m with a rate of 0.04 kg/ms. That means that we set the value for the Neumann boundary flux to that rate in between 1.75m and 2m. On the rest of the Neumann boundary we set the flux to zero, which means we define it as a no-flow boundary.
In addition, the DNAPL is injected at a point source at x = 0.502m and y = 3.02m with a rate of 0.1 kg/s.

We discretize the equations with a cell-centered finite volume TPFA scheme in space and an implicit Euler scheme in time. We use Newton's method to solve the system of nonlinear equations.
The grid is adaptively refined around the injection. The adaptive behaviour can be changed with input parameters in the `params.input` file.

# Implementation

## Folder layout and files

```
└── 2pinfiltration/
    ├── CMakeLists.txt          -> build system file
    ├── main.cc                 -> main program flow
    ├── params.input            -> runtime parameters
    ├── properties.hh           -> compile time configuration
    ├── problem.hh              -> boundary & initial conditions
    ├── spatialparams.hh        -> spatial parameter fields
    └── initialsolutioncc.txt   -> text file with initial solution
```
In the documentations behind the links provided in the following, you will find a detailed description how the the above mentioned set-up is realized in this example.
