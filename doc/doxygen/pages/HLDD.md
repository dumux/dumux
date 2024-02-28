# High-Level Design Documentation for DuMux

## Diagrams and Visual Aids

### Architecture Diagrams
- Diagrams representing the architecture and integration with DUNE and output possibilities.


@startuml
scale 0.8
top to bottom direction
package "Dune" {
package "Dune_ISTL"{
class ISTL_Solver
class ISTL_Matrices
class ISTL_Vectors
}
package DUNE_Grid{
class YASP
}

package "DuMux" {
package "Grid" {
interface GridManager{}
interface GridView{}
class GridGeometry{}
    }
    package "Variables" {
       class GridVariables{}
    }
    package "Scenario" {
        class Problem{}
        class TimeLoop{}
        class SolutionVector{}
    }
    package "Assembly" {
        class Assembler{}
        class LocalAssembler{}
        class LocalResidual{}
    }
    package "Solving" {
        class Solver{}
        interface LinearSolver{}
    }
    package "Output" {
        class IOField{}
        class VTKOutputModule
    }
}
}
LinearSolver <-- "Dune_ISTL"
GridManager <-- "DUNE_Grid"
GridManager --> GridView
GridView --> GridGeometry
GridGeometry --> Problem
GridGeometry --> GridVariables
GridGeometry --> Assembler
GridGeometry --> SolutionVector
TimeLoop --> Problem
LinearSolver --> Solver
LocalResidual --> LocalAssembler
LocalAssembler --> Assembler
GridVariables --> VTKOutputModule
GridVariables --> SolutionVector
IOField --> VTKOutputModule
Assembler --> Solver
SolutionVector --> Assembler
Problem --> GridVariables
Problem --> Assembler
GridVariables --> Assembler
}
@enduml

## Major Concepts and Components

In the following, a very brief overview over the major concepts and components of DuMux is provided. The focus is on the common ones that a user is exposed to for running a typical simulation. They are grouped into six categories: Grid, Variables, Assembly, Solving, Output, and Scenario. The interplay of the components is visualized in the diagram of the next section.

### Grid

#### GridManager
A GridManager abstracts and centralizes the creation, manipulation, and management of the computational grid to be used in the simulation. It serves as an interface to the underlying grid data structure coming in form of a DUNE grid implementation.

#### GridView
A GridView allows for read-only access to a certain part of a possibly hierarchical DUNE grid from which it is obtained. Most commonly employed is the LeafGridView, namely, a view on all elements of the grid without descendants in the hierarchy (which would be the leaves of a grid hierarchy).

#### GridGeometry
The GridGeometry constructs, from a GridView, all the geometrical and topological data necessary to evaluate the discrete equations. It is dependent on the selected spatial discretization method.


### Variables

#### GridVariables
GridVariables provide access to all variables needed to solve a particular discretized PDE, that is, the primary and secondary variables at geometric locations. These locations and also the type of variables depend on the chosen discretization scheme.


### Assembly

#### Assembler
The Assembler is responsible for calculating the global residual vector and the global Jacobian matrix. It relies on a discretization-specific local assembler engine. The Jacobian matrix is the partial derivative of the residual with respect to each entry of the solution vector. DuMux uses an element-wise assembly algorithm. For each element, a local assembler is instantiated.

#### LocalAssembler
The LocalAssembler is responsible for calculating the local residual vector and the local Jacobian matrix. The local residual vector is the vector of residuals for each element. The local Jacobian matrix is the partial derivative of the local residual with respect to each entry of the solution vector. The local assembler is discretization-specific.

#### LocalResidual
The LocalResidual is the implementation of a how a residual is evaluated on an element. It relies on the concept that each PDE has a storage, flux and source term. The storage and source terms are evaluated on the sub-control-volumes of an element and the flux term on the relevant faces of these sub-control-volumes.


### Solving

#### LinearSolver
The LinearSolver is a wrapper for a DUNE-ISTL preconditioned linear solver and used to solve the linear system of equations. It provides a common interface for different linear solvers.

#### Solver
The Solver manages the iterative refinement of solutions by assembling the Jacobian and residuals, solving the linearized equations, and applying the solution updates. Furthermore, it handles solution acceptance criteria.


### Output

#### IOFields
The IOFields class is responsible for managing the input and output fields. It provides methods and member functions for reading input files, initializing VTKOutputModules and managing the fields that are written to the output files.

#### VtkOutputModule
The VtkOutputModule is responsible for writing simulation results to VTK files for visualization. It can customize the output by adding variables to the output files. It generates one file per print-out step, possibly agglomerating several files from individual processes, and groups them into a PVD file containing time-step information.


### Scenario

#### Problem
A Problem in DuMux represents the conceptual framework where the scenario being simulated is characterized through the specification of initial and boundary conditions, as well as source terms.

#### TimeLoop
The TimeLoop manages temporal aspects of transient PDEs, handling parameters such as time-step size, current simulation time, and total simulation time. Stationary PDEs bypass the need for such temporal management.

#### SolutionVector
A SolutionVector is a container for the primary variables at each geometrical degree of freedom (dof). In particular, it holdsg NumEqVectors for each dof.

### User Interface
- Building your model is done using the Property system in DuMux
- Concretization of your model is done using the problem and params file
- Results can be exported to VTK