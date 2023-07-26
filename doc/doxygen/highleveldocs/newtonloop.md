@page newtonLoop Newton Loop

## Overview

A Newton algorithm is an iterative method used to solve non-linear equations, including non-linear partial differential equations (PDEs). It works by starting with an initial guess for the solution to the PDE, and then iteratively improving that guess according to their iterative changes until a solution is found.

Here's how it works in more detail:

1. Start with an initial guess for the solution to the PDE.
2. Approximate the partial derivatives for each DOF and assemble both the jacobian matrix and the residual [@ref assembler, @ref localassembler, @ref localresidual]
3. Solve linear system with the Jacobian and the residual to obtain a correction to the initial guess. [@ref linearsolver]
4. Correct the initial guess to obtain a new guess for the solution to the PDE.
5. Repeat steps 2-4 until the solution converges to a desired level of accuracy. (In the case of a linear PDE, only one iteration (steps 1-4) is required.

In addition, the NewtonSolver can recommend a time step for the next iteration of the time loop. In cases where more newton iterations were needed than expected, the recommended time step will be shorter, and in cases where fewer we needed, the recommended time step will increase.

##### Key functionalities

- solve()
   - Call of the assembleJacobianAndResidual() function of the @ref assembler
   - Call of the solve function of the @ref linearsolver
   - Update the solution using the newtonUpdate() function
   - Evaluate if convergence is reached using the newtonConverged() function
