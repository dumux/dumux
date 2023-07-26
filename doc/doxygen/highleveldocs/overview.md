# Overview

To solve your problem of interest, Dumux takes a PDE, transforms via a @ref Discretization to a set of algebraic equations. Afterwards the set of equations are solved.

On the highest level, Dumux only needs a few things to work:

1. First you must decide what kind of PDEs you will be solving, and with which discretization method is suited for this. (see \ref modelConcept \ref Discretization)
2. Once you specified your mathematical and numerical model concepts, you need to specify how the set of algebraic equation will be obtained. This is specified in a @ref assembler. The specification of the @ref assembler needs information about key properties:
   1. Geometry: All of the important locations and geometries needed for your discretization are stored in the @ref gridgeometry.
   2. Variables: everything related to primary and secondary variables is stored in the @ref gridvariables.
   3. LocalResidual: Your PDE of interest divided in Storage, Flux and Source Terms. See. @ref assembler
3. Information on how to solve the set of algebraic equations. Here either a See: @ref solver
4. Information about inital- and/or boundary-conditions. See: @ref problem
5. A container to store the solution of your problem. See: @ref solutionvector

From here on things are not mandatory, however in most cases you will encounter them:

6. If your problem is time dependent you need to specify a `TimeLoop`. See @ref timeloop
7. For output of your simulation you need to specify a `VTKOutputModule`

@mermaid{overview}


## Table of Contents

- @ref modelConcept
- @ref grid
   - @ref gridmanager
   - @ref gridview
   - @ref gridgeometry
- @ref gridvariables
   - @ref volumevariables
- @ref problem
- @ref solutionvector
- @ref assembler
- @ref solver
- @ref timeloop
- @ref vtkoutputmodule
   - @ref iofields
