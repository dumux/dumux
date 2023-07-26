# Solver

When it comes to solving the system, the PDE of interest must first classified as one of the following:

1. Non-linear PDEs: To solve a non-linear PDE, Dumux uses Newton's method to iteratively solve the PDE. The procedure is described below and in \ref newtonLoop.

2. Linear PDEs: To solve a linear PDE, a simplified version of Newton's method is developed as only one iteration is required. This is implemented as the `LinearPDESolver`, as discussed below.

Both cases additionally rely on a `LinearSolver` method. Both systems, linear or non-linear, produce a system of equations in the form of $\mathbf{J}\mathbf{x}=\mathbf{r}$, which need to be solved for each newton iteration (or at least once). Various `LinearSolver` methods are available for this step and must be specified (see also @ref linearsolver).

### Non-Linear PDE Solver (Newton Solver)

In order to iteratively solve non-linear PDEs, a newton's method algorithm is implemented and used in Dumux. This algorithm works by starting with an initial guess for the solution to the PDE, and then iteratively improving that guess until a solution is found that satisfies the given convergence criteria. For details on how this works, see \ref newtonLoop.

The `NewtonSolver` class is responsible for executing the Newton algorithm as described here.
It instantiates the assembly of the linear system and employs a `LinearSolver` to solve the system of equations (@ref linearsolver).
Depending on the complexity of the PDE and the convergence criteria, a number of newton iterations will be required before a solution is found.
For each iteration, the solution will be evaluated and updated, until convergence.

### Linear PDE Solver (Simplified Newton Solver)

`LinearPDESolver` implements a method used to solve linear PDEs.
This solver will do the same thing as the non-linear `Newton Solver` case, but when dealing with linear PDEs, only one iteration of the newton loop is necessary. As described in \ref newtonloop, an initial guess and a linear solver are required.

In addition, the jacobian matrix can be reused for further time steps as it will not change. This can be set using the function `reuseMatrix()`.
