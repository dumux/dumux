# Newton Loop

A Newton algorithm is an iterative method used to solve non-linear equations, including non-linear partial differential equations (PDEs). It works by starting with an initial guess for the solution to the PDE, and then iteratively improving that guess until a solution is found. The algorithm is implemented in {ref}`newtonSolver`.

Here's how it works in more detail:

1. Start with an initial guess for the solution to the PDE.
2. Use the PDE and the initial guess to form a linearized equation, called the Newton equation [{ref}`assembler`, {ref}`localAssembler`, {ref}`localResidual`]
3. Solve the Newton equation to obtain a correction to the initial guess.[{ref}`linearSolver`]
4. Add the correction to the initial guess to obtain a new guess for the solution to the PDE.
5. Repeat steps 2-4 until the solution converges to a desired level of accuracy.

## Overview
```{mermaid}
flowchart LR
    A(newtonSolver) -->|"assembleJacobianAndResidual()"| B(assembler)
    A -->|"solve()"| C(linearSolver)
    B -->|"assembleJacobianAndResidual()"| D(localAssembler)
    D -->|"evalFluxAndSourceResidual()"| E(localResidual)
    D -->|"evalStorageResidual()"| E
    click A "./newtonsolver.html"
    click B "./assembler.html"
    click C "./linearsolver.html"
    click D "./localassembler.html"
    click E "./localresidiual.html"
```
