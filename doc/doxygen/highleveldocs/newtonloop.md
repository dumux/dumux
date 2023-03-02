## NewtonLoop {#newtonloop}
<!-- @page newtonLoop Newton Loop -->

A Newton algorithm is an iterative method used to solve non-linear equations, including non-linear partial differential equations (PDEs). It works by starting with an initial guess for the solution to the PDE, and then iteratively improving that guess until a solution is found. The algortihm is implemented in @ref newtonSolver.

Here's how it works in more detail:

1. Start with an initial guess for the solution to the PDE.
2. Use the PDE and the initial guess to form a linearized equation, called the Newton equation [@ref assembler, @ref localassembler, @ref localresidual]
3. Solve the Newton equation to obtain a correction to the initial guess.[@ref linearsolver]
4. Add the correction to the initial guess to obtain a new guess for the solution to the PDE.
5. Repeat steps 2-4 until the solution converges to a desired level of accuracy.

## Overview

@mermaid{newtonloop}
