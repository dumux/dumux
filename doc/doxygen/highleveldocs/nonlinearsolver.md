@page nonlinearsolver NonLinearSolver

The nonLinearSolver class is used to execute the newton algorithm described in @subpage newtonloop "newtonloop". The assembly of the linear system is instantiated. Furthermore, a linear solver is used to solve the linear system of equations. In the newtonSolver class the evaluation of the solution is done.


### Key functionalities

1. solve()
   1. Call of the assembleJacobianAndResidual() function of the @ref assembler
   2. Call of the solve function of the @ref linearsolver
   3. Update the solution using the newtonUpdate() function
   4. Evaluate if convergence is reached using the newtonConverged() function

### Overview

@mermaid{newtonsolver}