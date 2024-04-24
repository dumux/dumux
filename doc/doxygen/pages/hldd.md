# High-Level Design Documentation for DuMux

## Diagrams and Visual Aids

### Architecture Diagrams
- Diagrams representing the architecture and integration with DUNE and output possibilities.

\pumlsvg{hldd,90}

## Major Concepts and Components

In the following, a brief overview over the major concepts and components of a DuMux simulation is provided. They are grouped hierarchically according to their code location and abstraction level, with the three main categories of User Scenario, DuMux and Dune. The interplay of the components is visualized in the diagram.

### User Scenario
On the top-most level, the user defines a scenario consisting of a problem and, often, associated spatial parameters. The physical model, the spatial discretization method and other compile-time choices are made by assigning respective properties. The `main` function is responsible for instantiating the necessary components, executing the simulation steps and writing the output files. Run-time parameters can be passed by a file like params.input.

#### Problem
A Problem in DuMux represents the conceptual framework where the scenario being simulated is characterized through the specification of initial and boundary conditions, as well as source terms.

#### SpatialParams
The SpatialParams define parameters which can be dependent on the location in space within the computational domain. For a Darcy-scale porous-medium simulation, these are typically porosity and permeability as well as parameters appearing in the constitutive relations, such as the capillary entry pressure.

#### Properties
In DuMux, properties are classes containing type definitions, values or methods which are selected at compile time. The selection includes, typically, the physical model, spatial discretization method, fluid system and grid manager. Choices are assembled in a file properties.hh.

#### main.cc
The main file includes the properties, by which, in particular, the problem is specified. In the `main` function, all necessary components are initiated and the simulation steps are performed.

#### params.input
The user can specify a file containing values for runtime parameters such as the number of grid cells or the initial time step size. The filename defaults to params.input.

### DuMux
The following lists the DuMux components that a user is exposed to for running a typical simulation. They are grouped into six categories: Grid, Variables, Assembly, Solving, Output, and Scenario.

#### Grid

##### GridManager
A GridManager abstracts and centralizes the creation, manipulation, and management of the computational grid to be used in the simulation. It serves as an interface to the underlying grid data structure coming in form of a DUNE grid implementation.

##### GridView
A GridView allows for read-only access to a certain part of a possibly hierarchical DUNE grid from which it is obtained. Most commonly employed is the LeafGridView, namely, a view on all elements of the grid without descendants in the hierarchy (which would be the leaves of a grid hierarchy).

##### GridGeometry
The GridGeometry constructs, from a GridView, all the geometrical and topological data necessary to evaluate the discrete equations. It is dependent on the selected spatial discretization method.


#### Variables

##### GridVariables
GridVariables provide access to all variables needed to solve a particular discretized PDE, that is, the primary and secondary variables at geometric locations. These locations and also the type of variables depend on the chosen discretization scheme.

##### PrimaryVariables
Vector type for storing the values of the independent variables at a geometric degree of freedom.

##### SolutionVector
A SolutionVector is a container for the primary variables at each geometrical degree of freedom (dof). In particular, it holdsg NumEqVectors for each dof.

#### Assembly

##### LocalResidual
The LocalResidual is the implementation of a how a residual is evaluated on an element. It relies on the concept that each PDE has a storage, flux and source term. The storage and source terms are evaluated on the sub-control-volumes of an element and the flux term on the relevant faces of these sub-control-volumes.

##### LocalAssembler
The LocalAssembler is responsible for calculating the local residual vector and the local Jacobian matrix. The local residual vector is the vector of residuals for each element. The local Jacobian matrix is the partial derivative of the local residual with respect to each entry of the solution vector. The local assembler is discretization-specific.

##### Assembler
The Assembler is responsible for calculating the global residual vector and the global Jacobian matrix. It relies on a discretization-specific local assembler engine. The Jacobian matrix is the partial derivative of the residual with respect to each entry of the solution vector. DuMux uses an element-wise assembly algorithm. For each element, a local assembler is instantiated.



#### Solving

##### LinearSolver
The LinearSolver is a wrapper for a DUNE-ISTL preconditioned linear solver and used to solve the linear system of equations. It provides a common interface for different linear solvers.

##### PDESolver
The PDESolver manages the iterative refinement of solutions by assembling the Jacobian and residuals, solving the linearized equations, and applying the solution updates. Furthermore, it handles solution acceptance criteria.

#### TimeLoop
The TimeLoop manages temporal aspects of transient PDEs, handling parameters such as time-step size, current simulation time, and total simulation time. Stationary PDEs bypass the need for such temporal management.


#### Output

##### IOFields
The IOFields class is responsible for managing the input and output fields. It provides methods and member functions for reading input files, initializing VTKOutputModules and managing the fields that are written to the output files.

##### VtkOutputModule
The VtkOutputModule is responsible for writing simulation results to VTK files for visualization. It can customize the output by adding variables to the output files. It generates one file per print-out step, possibly agglomerating several files from individual processes, and groups them into a PVD file containing time-step information.

### Dune

#### dune-common

##### ParameterTree
Implements a hierarchical structure of string parameters. Being accessible from practically any location within DuMux, it enables obtaining parameters passed by the user via a parameter file or the command line.

#### dune-grid

##### YaspGrid
The YaspGrid class is a structured, n-dimensional, parallel tensor product grid. It provides a distributed structured cube mesh and is designed to implement the DUNE grid interface for structured grids.

#### dune-istl

##### ISTL vectors
ISTL vector classes are designed to represent mathematical vector spaces. They support a recursive block structure, which is used to efficiently implement block preconditioners for hp-finite elements.

##### ISTL matrices
ISTL matrices classes are designed to represent linear maps between vector spaces. They also support a recursive block structure, which allows for efficient representation and computation.

##### ISTL solvers
ISTL solvers are designed to implement iterative solvers for linear systems in a generic manner. These solvers are subclasses of the abstract base class InverseOperator, which represents the inverse of an operator.
