# High-Level Documentation for DuMux

## General Description

### Product Perspective
DuMux, built on the DUNE framework, offers functionalities like finite volume discretizations and multi-domain framework for model coupling.

### Tools Used
Key technologies include modern C++17 and MPI for parallel computing.

### General Constraints and Assumptions
Discussion on the GPL-3.0 license terms and their implications.

## Design Details

### Main Design Features
- list the main classes and only one section (2-3 sentences) what this class is about

### Application and Technology Architecture
- Software architecture focusing on modularity and integration with DUNE framework.

### User Interface and User Experience
- Interaction through code for simulation setup and execution.
- Main interfaces
- vtu output?

## Diagrams and Visual Aids

### Architecture Diagrams
- Diagrams representing the architecture and integration with DUNE and output possibilities.

@mermaid{
flowchart TB
        subgraph DuMux
            subgraph Grid ["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Grid"]
            A(Grid) --> B(GridManager)

            B --> C(GridView)
            C --> D(GridGeometry)
            end

            subgraph Variables
            D --> F(GridVariables)
            end

            subgraph External
            D --> E(Problem)
            ETwo(TimeLoop) --> E
            D --> G(SolutionVector)
            style G fill:#f9f,stroke:#333,stroke-width:4px
            E --> F
            G -.-> F
            end

            subgraph Assembly ["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Assembly"]
            D --> H{Assembler}
            H <-.-> G
            E --> H

            ITwo(LocalAssembler) --> I(LocalResidual)
            H --> ITwo
            end

            subgraph Solver ["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Solver"]
            J(LinearSolver) --> K(Solver)

            K <-.-> H
            end

            subgraph Output
            L(IOField) -.-> M(VTKOutputModule)
            G -.-> M
            F -.-> M
            end
        end

        subgraph Dune
            Z((DUNE)) --> A
            Z --> J
        end
}

### Flowcharts
- Workflow of simulations and model setups.

## Conclusion


### Future Work
Planned future developments for DuMux.