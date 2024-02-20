# High-Level Design Documentation for DuMux

## General Description

### Product Perspective
DuMux, built on the DUNE framework, offers functionalities like finite volume discretizations, physics and multi-domain framework for model coupling.

### Tools Used
Key technologies include modern C++17 and MPI for parallel computing.

### General Constraints and Assumptions
Discussion on the GPL-3.0 license terms and their implications.

## Design Details

### Main Design Features
- list the main classes and only one section (2-3 sentences) what this class is about
#### Grid
#### GridManager
The concept of the GridManager in the DuMuX framework is to abstract and centralize the creation, manipulation, and management of various grid types used in numerical simulations. It serves as an interface between the user and the underlying grid data structures.
#### GridView
A GridView is a representation of the grid that allows for read-only access to certain parts of the Grid from which it is obtained. A leafGridView is a view on all elements of the grid without descendants in the hierarchy (which would be the leaves of a grid hierarchy) while a levelGridView is a view on all elements of a given level of a refinement hierarchy.
#### GridGeometry
The grid geometry constructs, from a GridView, all the geometrical and topological data necessary to evaluate the discrete equations.
#### GridVariables
GridVariables provide access to all variables needed to solve a particular (discrete) PDE, that is,
the primary and secondary variables at discrete locations.
These locations and also the type of variables can depend on the chosen discretization scheme.
#### Problem
The Problem class in DuMuX represents the conceptual framework where the scenario being simulated are characterized through the specification of initial and boundary conditions, as well as source terms.
#### TimeLoop
The TimeLoop class is instantiated to manage temporal aspects of transient PDEs, handling parameters such as time step size, current simulation time, and total simulation time, whereas stationary PDEs bypass the need for such temporal management.
#### SolutionVector
A SolutionVector class object is a container containing the primary variables for each degree of freedom (dof).In more detail, the SolutionVector object is a container class holding NumEqVector objects for each dof.
#### Assembler
The Assembler is responsible for calculating the global residual vector and the global Jacobian matrix. It relies on a discretization-specific local assembler engine. The Jacobian matrix is the partial derivative of the residual with respect to each entry of the solution vector. Dumux uses an element-wise assembly algorithm. For each element, a local assembler is instantiated.
#### LocalAssembler
The LocalAssembler is responsible for calculating the local residual vector and the local Jacobian matrix. The local residual vector is the vector of residuals for each element/scv. The local Jacobian matrix is the partial derivative of the local residual with respect to each entry of the solution vector. The local assembler is discretization-specific.
#### LocalResidual
The LocalResidual is the implementation of a how a residual is evaluated on a element/scv. It relies on the concept that each PDE has a storage,flux and source term. Here the actually PDE is implemented using this concept.
#### LinearSolver
The LinearSolver class is a wrapper for the actual linear solver used to solve the linear system of equations. It provides a common interface for different linear solvers. In the concrete linearSolver the implementation of the actual solving algorithm is written.
#### Solver
The Solver class manages the iterative refinement of solutions by assembling the Jacobian and residuals, solving the linearized equations, and applying the solution updates. Furthermore, it handles solution acceptance criteria.
#### IOField
#### VTKOutputModule
#### Coupling Manager

### Application and Technology Architecture
- Software architecture focusing on modularity and integration with DUNE framework.

### User Interface and User Experience
- Interaction through code for simulation setup and execution.
- Main interfaces
- vtu output?
- properties

## Diagrams and Visual Aids

### Architecture Diagrams
- Diagrams representing the architecture and integration with DUNE and output possibilities.

```plantuml
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
    package "External" {
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
}
```

### Flowcharts
- Workflow of simulations and model setups.

## Conclusion


### Future Work
Planned future developments for DuMux.
