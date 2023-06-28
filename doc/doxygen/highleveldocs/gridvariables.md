# GridVariables

`GridVariables` provide access to all variables needed to solve a particular (discrete) PDE, that is,
the primary and secondary variables at discrete locations.
These locations and also the type of variables can depend on the chosen discretization scheme.
For instance, `GridVariables` implementations for finite-volume schemes expose `GridVolumeVariables` and
a `GridFluxVariablesCache`. The former allow access to variables defined on sub-control volumes of the
finite-volume grid, while the latter exposes precomputed and cached variables that can be used to assemble
fluxes across sub-control volume faces.
Finally, note that for instationary problems, `GridVariables` expose two instances of `GridVolumeVariables`:
one for the current and another one for the volume variables of the previous time level.
### Key functionalities

- init(x): initalizes the grid variables for a given solution vector (primary variables) `x`
- update(x): updates the variables for a new solution `x`. Some variables, for instance the permeability in porous-medium flow problems, can be marked as solution-independent. In contrast to `init`, `update` only computes those variables that are solution-dependent, which can be significantly faster.


### Overview


@mermaid{gridvariables}

## GridVolumeVariables

As the name suggest, the `GridVolumeVariables` refer to variables that are located inside volumes on your grid.

In general, there are two options in Dumux:
1. Caching enabled
2. Caching disabled

If caching is enabled, a `GridVolumeVariables` instance (e.g., curGridVolVars) has a vector of `VolumeVariables`.
Since the smallest entity of a volume is a `Sub-control-volume` in Dumux, the length of the vector is the number of sub-control-volumes.
One can access entities of that vector with the respective sub-control-volume-index.

If caching is disabled, one can acces the `VolumeVariables` only via the `ElementVolumeVariables`.

Exemplary for the two implementations (e.g., caching enabled and disabled) in DuMuX, the implementation of `GridVolumeVariables` is sketched:
If caching is enabled, this is exactly how gridVolumeVariables are implemented. If caching is disabled, the global vector on the left does not exist.


<div align="center">
  <img src="GridVariables.svg"  width="50%">
</div>

### Key functionalities

If caching is enabled:

- update()
    - for every scv you will update the volumeVariables
- volVars()
    - gives you acces to volumeVariables stored at specific entry of the global gridVolumeVariables vector.

## ElementVolumeVariables

In short, `ElementVolumeVariables` are the view of a element towards all `VolumeVariables` located at sub-control-volumes in it's respective stencil.

In the case of enabled caching, the `ElementVolumeVariables` forward the entries belonging to that element from the globally stored vector described in `GridVolumeVariables`.

If caching is disabled, the `ElementVolumeVariables` do create a vector of `VolumeVariables` for each element on the fly.

In both cases the `ElementVolumeVariables` will be forwarded to the `LocalAssembler`. Since the `LocalAssembler` only needs information of the sub-control-volumes in the respective stencil.

### Key functionalities

- bind()
  - for caching enabled, precomputes dirichlet boundary values in the stencil
  - for caching disabled, computes all @ref volumevariables in the stencil.

## VolumeVariables

`VolumeVariables` is a class that gives access and lets you store variables are actually needed for the computation. Here, one could implement a function that calculates the density of your fluid of interest from `PrimaryVariables`.
However, in DuMuX there exists another layer. For instance, access to variables that are connected to your fluid of interest can be accessed and stored via the fluidstate.
If you keep using the example of variables that are connected to your fluid of interest, equations to calculate the variables are per default defined in the  fluidsystem.
The same logic applies to other scenarios, for instance variables that are connected to the solid (i.e., porosity).
Per default, the VolumeVariables class acts gives you access to the next layer. However, it is not mandatory to implement a `FluidState` or `Fluidsystem`.

For implementation details it is reffered to the modules-documentation:

- @ref FluidStates
- @ref FluidSystems
- @ref SolidStates
- @ref SolidSystems

### Key functionalites

- completeFluidState()
  - sets the properties of the @ref fluidstate, that are primaryVariables
  - calculates the rest of the variables concerning the fluid using the functions defined in the @ref fluidsystem