@page linearpdesolver LinearPDESolver

As the name suggest, this is a solver that solve linear PDE's.
It does the same thing as the @ref nonlinearsolver, but for linear PDE's only one iteration is necessary.
The assembly of the linear system is instantiated. Furthermore, a linear solver is used to solve the linear system of equations.

### Key functionalities

1. solve()
   1. Call of the assembleJacobianAndResidual() function of the @ref assembler
   2. Call of the solve function of the @ref linearsolver

### Overview

@mermaid{linearpdesolver}