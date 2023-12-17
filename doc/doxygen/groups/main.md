# Doxygen top level groups

@defgroup Models Mathematical models
@brief The mathematical models implemented in DuMu<sup>x</sup>

@defgroup Discretization Discretization schemes
@brief The discretization schemes available in DuMu<sup>x</sup>

@defgroup Properties Properties and the property systems

@defgroup Parameter Parameters and runtime configuration

@defgroup Material Constitutive models
@brief Constitutive relations for fluids, solids, fluid-matrix interactions, and more
@details Constitutive relations formalize the functional dependence among physical variables, thereby providing the missing equations to close systems. Some constitutive relations, in particular if they are ad-hoc (e.g. Hooke's law) or empirical
and not derived from basic principles, are also called material laws.
Implementations of the functional relationships available in DuMu<sup>x</sup> are collected in the folder `dumux/material`.

Most implementations follow one of the following paradigms

* free functions (in nested namespaces),
* stateless classes with static member functions where all parameters are function arguments,
* statefull classes that contain parameters that need to be initialized.

Stateless, purely static classes have the advantage that only the class type needs to be known to
the class using them. Statefull classes need to be instantiated
and the instantiated object is passed to the consumer class.

To view the description of individual constitutive models and classes, click on one the modules below.

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
