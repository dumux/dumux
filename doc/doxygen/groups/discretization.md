# Doxygen subgroups of the discretization group

@defgroup CCDiscretization Cell-centered Finite Volume Methods
@brief Finite volume schemes with degrees of freedom located at grid cell centers.
@ingroup Discretization

@defgroup CCTpfaDiscretization Two-point flux approximation (Tpfa)
@brief A cell-centered finite volume scheme with two-point flux approximation.
@ingroup CCDiscretization

@defgroup CCMpfaDiscretization Multi-point flux approximation (Mpfa)
@brief A cell-centered finite volume scheme with multi-point flux approximation.
@ingroup CCDiscretization

@defgroup CVFEDiscretization Control-volume Finite Element Methods
@brief Control-volume finite element schemes (e.g. box method)
Control-volume finite element schemes are based on finite element basis functions for interpolation but define control volumes to construct a finite volume scheme. They can be interpreted both as finite volume or as (Petrov-Galerkin) finite element scheme.
@ingroup Discretization

@defgroup DiamondDiscretization Diamond discretization scheme
@brief Face-centered finite-volume scheme based on non-conforming finite-element spaces.
@ingroup CVFEDiscretization

@defgroup BoxDiscretization Box FV scheme
@brief The box method is a collocated finite volume scheme with control volumes centered at grid nodes.
@ingroup CVFEDiscretization

@defgroup PQ1BubbleDiscretization PQ1 bubble scheme
@brief Control-volume finite element scheme based on P1/Q1 basis function with enrichment by a bubble function
@ingroup CVFEDiscretization

@defgroup FaceCenteredStaggeredDiscretization Staggered Grid Finite Volume Method
@brief Discretization for the momentum balance of the Navier-Stokes equations. Can be used to build a marker-and-cell scheme (MAC) together with Tpfa for the discretization of the mass balance equation.
@ingroup Discretization

@defgroup StaggeredDiscretization Staggered FV scheme
@brief A staggered finite volume scheme with degrees of freedom at cell-centers and facets. In this implementation, momentum control volumes do not explicitly exist, but the implementation uses workarounds.
@note This is an outdated implementation of the MAC scheme and will be replaced. See @ref FaceCenteredStaggeredDiscretization
@ingroup Discretization

@defgroup FEMDiscretization Finite Element Methods
@brief The finite element method
@ingroup Discretization

@defgroup PoreNetworkDiscretization Pore-network Models
@brief The pore-network model discretization.
@ingroup Discretization
