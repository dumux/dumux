@page linearsolver LinearSolver

## LinearSolver

The `LinearSolver` solves a linear-system of equations of the form $\mathbf{J}\mathbf{x}=\mathbf{r}$,
where $\mathbf{J}$ denotes the Jacobian while $\mathbf{r}$ represents the residual.

The remaining $\mathbf{x}$ is then used as the correction to the previous solution. $mathbf{x_{new}} = mathbf{x_{old}} - mathbf{x}$

The used linear solvers are implemented in DUNE. Access to these linear solvers are developed in Dune-ISTL (Iterative Solver Template Library)
<a href="https://www.dune-project.org/modules/dune-istl/">Dune-ISTL (Iterative Solver Template Library)</a>
