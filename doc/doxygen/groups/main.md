# Doxygen top level groups

@defgroup Models Mathematical models
@brief The mathematical models implemented in DuMu<sup>x</sup>

@defgroup Discretization Discretization schemes
@brief The discretization schemes available in DuMu<sup>x</sup>

@defgroup Properties Properties and the property systems

@defgroup Parameter Parameters and runtime configuration

@defgroup Material Constitutive modeling
@brief Constitutive model framework, material models, fluids, solids.

@defgroup MultiDomain Multidomain framework
@brief Coupling of several regular DuMu<sup>x</sup> problems
@details The multi domain module allows coupling regular DuMu<sup>x</sup> problems.
Several coupling modes are currently available.

@defgroup Geometry Geometry
@brief Algorithms for geometry computations (intersections, distances, ...).

@defgroup InputOutput Input Output
@brief Input and output of data and grids

@defgroup AssemblyAndSolvers Assembly and Solvers
@brief Assembling matrices and vectors, solvers for linear and nonlinear equations

<!-- AssemblyAndSolvers subgroups begin -->

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

<!-- AssemblyAndSolvers subgroups end -->

@defgroup Typetraits Typetraits
@brief Basic Type traits in DuMu<sup>x</sup>

@defgroup Core Common
@brief Common functionality

@defgroup Flux Numerical flux approximations
@brief Everything flux related in DuMu<sup>x</sup>

@defgroup SpatialParameters Spatial parameters
@brief Spatial parameters.
@details Spatial parameters.

@defgroup Adaptive Tools for adaptive grids
@brief Tools for simulations using adaptive grids

@defgroup Experimental Experimental
@brief Experimental features
