# Solver

Depending on the PDE of interest you have either one of the two kinds of PDE's:

1. linear PDE: To solve a linear PDE a `LinearPDESolver` is implemented.
2. non-linear PDE: To solve a non-linear PDE, Dumux uses Newton's method to solve the PDE. The procedure is described below. To solve a non-linear PDE in Dumux you need to create an `NewtonSolver` instance.
3. Both methods additionaly rely on a `LinearSolver`, since both methods produce a system of equations in the form Ax=b. To solve that system you have to specify a `LinearSolver`.

## NewtonSolver
The Newton algorithm is an iterative method used to solve non-linear equations, including non-linear partial differential equations (PDEs). It works by starting with an initial guess for the solution to the PDE, and then iteratively improving that guess until a solution is found. The algorithm is implemented in `NewtonSolver`.

Here's how it works in more detail:

1. Start with an initial guess for the solution to the PDE.
2. Use the PDE and the initial guess to form a linearized equation, called the Newton equation
3. Solve the Newton equation to obtain a correction to the initial guess.
4. Add the correction to the initial guess to obtain a new guess for the solution to the PDE.
5. Repeat steps 2-4 until the solution converges to a desired level of accuracy.

The `NewtonSolver` class is responsible for executing the Newton algorithm as described earlier. It instantiates the assembly of the linear system and employs a `LinearSolver` to solve the system of equations. Furthermore, evaluating and updating you solution is done inside the `NewtonSolver`.

### Key functionalities

- solve()
   - Call of the assembleJacobianAndResidual() function of the `Assembler`
   - Call of the solve function of the `LinearSolver`
   - Update the solution using the newtonUpdate() function
   - Evaluate if convergence is reached using the newtonConverged() function

## LinearPDESolver

As the name suggest, this is a solver that solve linear PDE's.
It does the same thing as the `NewtonSolver`, but for linear PDE's only one iteration is necessary.

### Key functionalities

1. solve()
   1. Call of the assembleJacobianAndResidual() function of the @ref assembler
   2. Call of the solve function of the @ref linearsolver
## LinearSolver

The `LinearSolver` solves a linear-system of equations in  the form Ax=b.
In the case of Dumux it is Jx=r, with the Jacobian J and the residual r.
The used linear solvers are implemented in DUNE.


### Key functionalities

- solve()
  - solves a system of equation of the Form Jx=r
