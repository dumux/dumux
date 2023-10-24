# Doxygen subgroups of the multidomain group

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
