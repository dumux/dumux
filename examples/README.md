# Examples

To get started with DuMu<sup>x</sup>, we recommend the following documented examples. Each example folder contains a ready-to-use DuMu<sup>x</sup> simulation example.
The description in each folder (best viewed online by following the link) explains each line of the code example.

### [:open_file_folder: Example 1: One-phase flow and tracer transport](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples/1ptracer)

In this example, we simulate tracer transport through a confined aquifer with a randomly distributed permeability field.
We first solve the pressure field, compute the steady state flow field, and then solve the tracer transport equation.
You learn how to

* generate a randomly distributed permeability field
* solve a one-phase flow in porous media problem
* compute the flow field from a pressure solution to pass to a tracer problem
* sequentially solve two types of problems after each other

### [:open_file_folder: Example 2: Two-phase flow with infiltration and adaptive grid](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples/2pinfiltration)

In this example we model a soil contamination problem where DNAPL infiltrates a water-saturated porous medium (two-phase flow).
The initial distribution of DNAPL is read in from a txt-file.
The grid is adapitvely refined where DNAPL enters the domain, around the plume, and around an injection well.
You learn how to

* solve a two-phase flow in porous media problem with two immiscible phases
* set boundary conditions and a simple injection well
* implement a problem with heterogenous material parameters
* use adaptive grid refinement around the saturation front

### [:open_file_folder: Example 3: Shallow water model](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples/shallowwaterfriction)

The shallow water flow model is applied to simulate steady subcritical flow in a river including a bottom friction model.
You learn how to

* solve a shallow water flow problem including bottom friction
* computate and output (VTK) an analytical reference solution

### [:open_file_folder: Example 4: Freeflow channel](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples/freeflowchannel)
In this example, we simulate a free flow between two plates in two dimensions.
You learn how to

* solve a free flow problem
* set outflow boundary conditions in the free-flow context
