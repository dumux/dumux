# Basic concepts

## Properties / Property System
In DuMuX, the property system provides a flexible way to configure simulations at compile time. Properties are structs that define types which determine how different parts of the framework work together. The system ensures a consistent class hierarchy by allowing users to inherit from predefined models and customize specific properties as needed.
Users typically collect their property customizations in a `properties.hh` file specific to their simulation setup. For more detailed and technical information about the property system, see @ref Properties.

## Problem
A Problem in DuMux represents the conceptual framework where the scenario being simulated is characterized through the specification of initial and boundary conditions, as well as source terms. It is recommended to save all specifications of the initial and boundary conditions in one file e.g. `problem.hh`.

## SpatialParams
The SpatialParams define parameters which can be dependent on the location in space within the computational domain. For a Darcy-scale porous-medium simulation, these are typically porosity and permeability as well as parameters appearing in the constitutive relations, such as the capillary entry pressure.The specified spatial parameters are usually stored in a `spatialparams.hh` file. More information can be found @ref SpatialParameters.

## Grid
### GridManager
A GridManager abstracts and centralizes the creation, manipulation, and management of the computational grid to be used in the simulation. It serves as an interface to the underlying grid data structure coming in form of a DUNE grid implementation.

### GridView
A GridView allows for read-only access to a certain part of a possibly hierarchical DUNE grid from which it is obtained. Most commonly employed is the LeafGridView, namely, a view on all elements of the grid without descendants in the hierarchy (which would be the leaves of a grid hierarchy).


### GridGeometry
The GridGeometry constructs, from a GridView, all the geometrical and topological data necessary to evaluate the discrete equations. It is dependent on the selected spatial discretization method.

## Assembly

### LocalResidual
The LocalResidual is the implementation of a how a residual is evaluated on an element. It relies on the concept that each PDE has a storage, flux and source term. The storage and source terms are evaluated on the sub-control-volumes of an element and the flux term on the relevant faces of these sub-control-volumes.

### LocalAssembler
The LocalAssembler is responsible for calculating the local residual vector and the local Jacobian matrix. The local residual vector is the vector of residuals for each element. The local Jacobian matrix is the partial derivative of the local residual with respect to each entry of the solution vector. The local assembler is discretization-specific.

### Assembler
The Assembler is responsible for calculating the global residual vector and the global Jacobian matrix. It relies on a discretization-specific local assembler engine. The Jacobian matrix is the partial derivative of the residual with respect to each entry of the solution vector. DuMux uses an element-wise assembly algorithm. For each element, a local assembler is instantiated.

## Solving

### LinearSolver
The LinearSolver is a wrapper for a DUNE-ISTL preconditioned linear solver and used to solve the linear system of equations. It provides a common interface for different linear solvers.

### PDESolver
The PDESolver manages the iterative refinement of solutions by assembling the Jacobian and residuals, solving the linearized equations, and applying the solution updates. Furthermore, it handles solution acceptance criteria.
