@addtogroup MultiDomain

![Coupled models realized with the DuMux multidomain framework](multidomain.svg){html: width=100%}

DuMu<sup>X</sup> can couple problems posed on different domains.
The domains can touch or overlap, model different physics, have different dimensions, different grids, or different discretization methods.
The full system Jacobian is approximated by numeric differentiation (avoids hand-coding the Jacobian)
which allows building monolithic solvers for complex nonlinear coupled problems.
The framework enabling all these coupled simulations is called the MultiDomain framework.

## Design Goals

* Reuse existing DuMux models in the coupled setting
* Extensible to more than two subdomains
* Common structure for different coupling modes

## Monolithic assembly and solvers

![Structure of monolithic Jacobian matrix](mdstructure.png){html: width=50%}

The multidomain framework uses block-structured matrices to assemble a single system matrix for the linear system
arising from the coupled problem (for example, in each Newton iteration of a nonlinear solver).

The matrix is structured such that each of the diagonal blocks correspond to one of the subproblems.
The off-diagonal blocks correspond to the coupling derivatives that represent interactions between the subproblems.
The block-structured matrix can directly be used with compatible linear solvers.
The block structure is also useful to construct partitioned solvers.

## The Subproblems

Any regular ("single domain") DuMux problem can be turned into a subproblem of a coupled simulation.
The subproblems can model arbitrary physics, can have different dimensions, or be discretized with different methods.
The only requirement is that the subproblems are compatible with the coupling manager.

In addition to a regular DuMux problem, the subproblems gets a pointer to the coupling manager which
is used to transfer data between the subproblems and obtain data from other subproblems.

## The Coupling Manager

The coupling manager is the central object in the multidomain framework,
which is responsible for the coupling of the subproblems.
Coupling managers implement the interface of the Dumux::CouplingManager class.

The coupling manager is responsible for the following tasks:

* Transfer data from one subproblem to another
    - e.g. give me the soil pressure (from the well domain)
    - e.g. give me the radius of the embedded well (from the soil domain)

* Compute the coupling stencil
* Compute the coupling residual (numerical differentiation)


@defgroup BoundaryCoupling Boundary coupling mode
@brief Couples problems of different or equal dimension that touch at the domain boundary. Examples are equal-dimension multi-physics problems like Darcy-Stokes coupling or PNM (pore network model)-Darcy coupling.
@ingroup MultiDomain

@defgroup DarcyDarcyCoupling Darcy-Darcy domain coupling
@brief Couples domains with equal-dimension multi-physics problems in a Darcy-Darcy coupling.
@ingroup BoundaryCoupling

@defgroup FreeFlowPoreNetworkCoupling Free flow-Pore network domain coupling
@brief Couples domains with equal-dimension multi-physics problems in a Free flow-Pore network coupling.
@ingroup BoundaryCoupling

@defgroup FreeFlowPorousMediumCoupling Free flow-Porous medium domain coupling
@brief Couples domains with equal-dimension multi-physics problems in a Free flow-Porous medium coupling.
@ingroup BoundaryCoupling

@defgroup StokesDarcyCoupling Stokes-Darcy domain coupling
@brief Couples domains with equal-dimension multi-physics problems in a Stokes-Darcy coupling.
@ingroup BoundaryCoupling

@defgroup DualNetworkCoupling Coupling for dual network approach for pore network models
@brief Coupling for dual network approach for pore network models
@ingroup MultiDomain

@defgroup EmbeddedCoupling Embedded mixed-dimension coupling mode
@brief Couples problems of different dimensions where one or more lower-dimensional problems (lowdim) are embedded in a higher-dimensional domain (bulk). Examples are embedded one-dimensional networks for the simulation of blood tissue perfusion, or root-soil interaction, and embedded fracture models.
@ingroup MultiDomain

@defgroup FacetCoupling Conforming mixed-dimension facet coupling mode
@brief Couples problems of different dimensions where one or more lower-dimensional problems (lowdim) live on the facets of the higher-dimensional domain (bulk). Examples are discrete facet conforming fracture models and problems with physics on a domain surface.
@ingroup MultiDomain
