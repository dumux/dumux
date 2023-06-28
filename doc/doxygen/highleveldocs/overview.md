# Overview

To solve your problem of interest, Dumux takes a PDE, transforms via a @ref Discretization to a set of algebraic equations. Afterwards the set of equations are solved.

On the highest level, Dumux only needs a few things to work:

1. Specify what type of solution strategy you want to apply. See: @ref solver
2. Once you specified your solver, you need to specify how the set of algebraic equation will be obtained. This is specified in a @ref assembler. The specification of the @ref assembler needs information about key properties:
   1. Geometry: everything related to Geometry is stored in the @ref gridgeometry.
   2. Variables: everything related to variables is stored in the @ref gridvariables.
   3. LocalResidual: basically your PDE of interest divided in Storage, Flux and Source Terms. See. @ref assembler
3. Information on how to solve the set of algebraic equations. Therefore, you must specify a `LinearSolver`. See: @ref solver
4. Information about inital- and/or boundary-conditions. See: @ref problem
5. A container to store the solution of you problem. See: @ref solutionvector

From here on things are not mandatory, however in most cases you will encounter them:

6. If your problem is time dependet you need to specify a `TimeLoop`. See @ref timeloop
7. For output of your simulation you need to specify `IOFields` and a `VTKOutputModule`

@mermaid{overview}


## Table of Contents

- @ref grid
- @ref gridmanager
- @ref gridview
- @ref gridgeometry
- @ref gridvariables
- @ref assembler
- @ref solver
- @ref timeloop
- @ref problem
- @ref solutionvector
- @ref iofields
- @ref vtkoutputmodule