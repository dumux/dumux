# Doxygen subgroups of the core group

<!-- Properties -->

@defgroup Properties Properties and the property systems
@ingroup Core

<!-- Parameters -->

@defgroup Parameters Parameters and runtime configuration
@ingroup Core

<!-- Discretization -->

@defgroup Discretization Discretization schemes
@brief The discretization schemes available in DuMu<sup>x</sup>
@ingroup Core

<!-- Discretization subgroups -->

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

<!-- Flux -->

@defgroup Flux Numerical flux approximations
@brief Everything flux related in DuMu<sup>x</sup>
@ingroup Core

<!-- Flux subgroups -->

@defgroup BoxFlux Flux related to the box scheme
@brief Flux related to the box scheme
@ingroup Flux

@defgroup CVFEFlux Flux related to the CVFE scheme
@brief Flux related to control-volume finite element schemes
@ingroup Flux

@defgroup CCFlux Flux related to the cell-centered schemes
@brief Flux related to the cell-centered schemes
@ingroup Flux

@defgroup CCTpfaFlux Flux related to the cell-centered two-point flux approximation schemes
@brief Flux related to the cell-centered two-point flux approximation schemes
@ingroup Flux

@defgroup CCMpfaFlux Flux related to the cell-centered multi-point flux approximation schemes
@brief Flux related to the cell-centered multi-point flux approximation schemes
@ingroup Flux

@defgroup FaceCenteredDiamondFlux Flux related to the face-centered diamond scheme
@brief Flux related to the face-centered diamond scheme
@ingroup Flux

@defgroup PoreNetworkFlux Flux related to the pore network models
@brief Flux related to the pore newtwork models
@ingroup Flux

@defgroup StaggeredFlux Flux related to the staggered scheme
@brief Flux related to the staggered scheme
@ingroup Flux

@defgroup ShallowWaterFlux Flux related to the shallow water model
@brief Flux related to the shallow water model
@ingroup Flux

<!-- SpatialParameters -->

@defgroup SpatialParameters Spatial parameters
@brief Spatial parameters.
@details Spatial parameters.
@ingroup Core

<!-- Adaptive -->

@defgroup Adaptive Tools for adaptive grids
@brief Tools for simulations using adaptive grids
@ingroup Core

<!-- AssemblyAndSolvers -->

@defgroup AssemblyAndSolvers Assembly and Solvers
@brief Assembling matrices and vectors, solvers for linear and nonlinear equations
@ingroup Core

<!-- AssemblyAndSolvers subgroups -->

@defgroup Assembly Assembly
@brief Assembly of linear systems (Jacobian and residual)
@ingroup AssemblyAndSolvers

@defgroup Linear Linear
@brief Linear solvers and helpers
@ingroup AssemblyAndSolvers

@defgroup Nonlinear Nonlinear
@brief Nonlinear solvers: Newton method
@ingroup AssemblyAndSolvers

@defgroup Parallel Parallel
@brief Files for communication of parallel solvers
@ingroup AssemblyAndSolvers

<!-- Geometry -->

@defgroup Geometry Geometry
@brief Algorithms for geometry computations (intersections, distances, ...).
@ingroup Core

<!-- InputOutput -->

@defgroup InputOutput Input Output
@brief Input and output of data and grids
@ingroup Core

<!-- Typetraits -->

@defgroup Typetraits Typetraits
@brief Basic Type traits in DuMu<sup>x</sup>
@ingroup Core
