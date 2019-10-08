Examples
===============
To get started with DuMu<sup>x</sup>, we recommend the following documented examples. Each example folder contains a ready-to-use DuMu<sup>x</sup> simulation example.
The description in each folder (best viewed online by following the link) explains each line of the code. 

 * [Example 1: One phase flow with a tracer model:](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples/1ptracer)

    In this example we simulate tracer transport through a confined aquifer with a randomly distributed permeability field. We first solve the pressure field, compute the steady state flow field,
    and then solve the tracer transport equation.

 * [Example 2: Two-phase flow with infiltration and adaptive grid:](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples/2pinfiltration)

    In this example we model a soil contamination problem where DNAPL infiltrates a water-saturated porous medium (two-phase flow).
    The initial distribution of DNAPL is read it from a txt-file.
    The grid is adapitvely refined where DNAPL enters the domain, around the plume, and around an injection well.

 * [Example 3: Shallow water model:](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples/shallowwaterfriction)

    The shallow water flow model is applied to simulate steady subcritical flow in a river including a bottom friction model.
