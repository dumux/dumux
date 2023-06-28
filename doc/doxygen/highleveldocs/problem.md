# Problem

Within the `Problem` class, the initial conditions, boundary conditions and volumetric sink/sources are defined. If spatial parameters are necessary, the `Problem` class also stores a pointer to the respective @subpage spatialparams "SpatialParams" object. Any custom functions, such as pre- or post-timestep operations could also be included here. All instances of a `Problem` class inherit from the [base Problem class](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/common/fvproblem.hh) or from [this base Problem class](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/common/fvproblemwithspatialparams.hh) if spatial parameters are also required..

### Key functionalities

* name():
    - Return the problem name.
* setName():
    - Set the problem name.
* boundaryTypes(element, scv/scvf):
    - Specifies which kind of boundary condition should be used for which equation for a given `element/scv/scvf`.
* boundaryTypesAtPos(pos):
    - Specifies which kind of boundary condition should be used for which equation for a given position `pos`.
* dirichlet(element, scv/scvf):
    - Evaluate the boundary conditions for a dirichlet `element/scv/scvf`.
* dirichletAtPos(pos):
    - Evaluate the boundary conditions for a dirichlet position `pos`.
* neumann(element, scvf, ...):
    - Evaluate the boundary conditions for a neumann boundary segment at `element/scvf`.
* neumannAtPos(pos):
    - Evaluate the boundary conditions for a neumann boundary segment at `pos`.
* source(element, scv, ...):
    - Evaluate the source term for all phases within a given sub-control-volume `scv`.
* sourceAtPos(pos):
    - Evaluate the source term for all phases at a given position `pos`.
* pointSource(element, scv, ..):
    - Evaluate the point sources for all phases within a given sub-control-volume `scv`.
* pointSourceAtPos():
    - Evaluate the point sources for all phases at a given position `pos`.
* applyInitialSolution(sol):
    - Applies the initial solution for all degrees of freedom of the grid given the solution vector `sol`.
* intial(entity):
    - Evaluate the initial value for an `entity`. Entities can be `elements`, `scvs`, `scvfs` and other entities describing the grid.
* initialAtPos(pos):
    - Evaluate the initial value for a given position `pos`.
* gridGeometry():
    - Return the finite volume grid geometry belonging to this `problem`.
* paramGroup():
    - Return the parameter group from which to retrieve runtime parameters.
