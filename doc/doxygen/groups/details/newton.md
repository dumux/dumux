@defgroup Newton Newton solver
@brief The Newton solver and the related parameters. This describes the reference implementation @ref Dumux::NewtonSolver.
@brief The Newton solver and the related parameters. This describes the reference implementation @ref Dumux::NewtonSolver.
@ingroup Nonlinear
@details

# Newton's method

The following describes Newton's method and the reference implementation @ref Dumux::NewtonSolver of a multi-dimensional Newton solver to solve a non-linear system of equations in DuMu<sup>x</sup>.
For solving nonlinear equations with only one unknown, DuMu<sup>x</sup> also provides @ref Dumux::findScalarRootNewton and @ref Dumux::findScalarRootBrent.

We assume that the equations (in DuMu<sup>x</sup> usually the discretized partial differential equations) are formulated as a residual equation:

```math
F(\textbf{u}) = \textbf{r} = \textbf{0}
```

with the residual operator $F : \mathbb{R}^n \to \mathbb{R}$, $n$ is the number of equations, $\textbf{u} \in \mathbb{R}^n$ is the vector of unknowns,
and $\textbf{r} \in \mathbb{R}^n$ is the residual vector.

The goal of Newton's method is to find a vector $\textbf{u}$ such that the residual equation is fulfilled (at least to a good enough approximation $F(\textbf{u}) \approx \textbf{0}$ as we will discuss below).

Newton's method is an iterative algorithm: it generates a sequence of approximations, $\textbf{u}_k$, that converges (under some circumstances) to the vector $\textbf{u}_\star$, that satisfies $F(\textbf{u}_\star) = 0$. We start with an initial guess $\textbf{u}_0$ and
calculate the initial residual $\textbf{r}_0 = \textbf{F}(\textbf{u}_0)$.
Then, we calculate the derivative of the residual operator with respect to the unknowns, the Jacobian matrix, $J : \mathbb{R}^n \to \mathbb{R}^n$
evaluated at $\textbf{u}_k$ (in the first step $k=0$):

```math
J(\textbf{u}_k) = \frac{\partial F(\textbf{u})}{\partial \textbf{u}} \vert_{\textbf{u} = \textbf{u}_k}
```

@note Computing the residual operator $F$ and the Jacobian matrix $J$ is delegated to the assembler implementation (see e.g. @ref Dumux::FVAssembler).
Typically with DuMu<sup>x</sup>, the Jacobian matrix is approximated by finite differences but the implementation of Newton's method in @ref Dumux::NewtonSolver is independent of the way the Jacobian is calculated. The C++ type of the Assembler is a template parameter of @ref Dumux::NewtonSolver and and instance of the assembler must be passed to the Newton solver in its constructor @ref Dumux::NewtonSolver::NewtonSolver.

This allows to compute a tangent direction of residual surface in n-dimensional space as

```math
\Delta \textbf{u} = \textbf{u}_{k} - \textbf{u}_{k+1} = J(\textbf{u}_k)^{-1} \cdot \textbf{r}_k
```

@note Evaluating $\Delta \textbf{u}$ is equivalent to solving the linear system
```math
J(\textbf{u}_k) \Delta \textbf{u} = \textbf{r}_k
```
This task is delegated to a linear solver (see e.g. @ref Dumux::LinearSolver)
which is passed the Jacobian matrix and the residual vector computed by the assembler earlier.
The C++ type of the LinearSolver is a template parameter of @ref Dumux::NewtonSolver and and instance of the linear solver must be passed to the Newton solver in its constructor @ref Dumux::NewtonSolver::NewtonSolver.

With the tangent direction $\Delta \textbf{u}$ (also called the "update" vector), we update our
current approximation $\textbf{u}_k$ to get a new approximation $\textbf{u}_{k+1}$:

```math
\textbf{u}_{k+1} = \textbf{u}_k - \Delta \textbf{u}
```

@note The update can be customized by a user-defined update strategy or by choosing from the implemented strategies (see below for more details).

This is one Newton step (or iteration). We repeat the steps of calculating the residual
and the Jacobian matrix, solve the linear system, and update the approximation until we reach a good enough solution (see below for the available termination criteria).

@note A good initial guess is crucial to get fast convergence of Newton's method.
For time-dependent problems, the initial guess $\textbf{u}_0$ is usually the solution of the previous time step.
The initial guess is passed to the Newton solver in the method @ref Dumux::NewtonSolver::apply and @ref Dumux::NewtonSolver::solve.

@note When choosing a zero initial guess you must be aware that the assembler (based on finite difference approximations)
may have difficulties to find a good epsilon to compute the finite difference approximation. In this case, you may need to
help the assembler by providing an estimate of the problem-typical magnitude of the unknowns. See @ref Dumux::FVAssembler
and the parameters `Assembly.NumericDifference.PriVarMagnitude` and `Assembly.NumericDifference.BaseEpsilon`.

The @ref Dumux::NewtonSolver reference implementation of Newton's method has several features discussed below that can be controlled by user-set parameters (e.g. via a parameter input file or the command line, see @ref runtime-parameters).

## Update strategy

Here we describe strategies that differ from the default update strategy, that is:

```math
\textbf{u}_{k+1} = \textbf{u}_k - \Delta \textbf{u}
```
The different strategies can be enabled via runtime parameters as described below.


### Line search

The idea of line search is to find an optimal scaling factor $\alpha \in \mathbb{R}$ for the update vector $\Delta \textbf{u}$
```math
\textbf{u}_{k+1} = \textbf{u}_k - \alpha\Delta \textbf{u}
```

such that the residual norm evaluated with the new iterate $\textbf{u}_{k+1}$ that is $||r_{k+1}||$ is minimized.
That means using line search the following optimization problem is solved: find $\alpha$ such that

```math
\alpha = \arg \min_{\alpha} ||F(\textbf{u}_{k+1} - \alpha \Delta \textbf{u})||_2^2
```

The default implementation of @ref Dumux::NewtonSolver uses a simple backtracking line search algorithm to approximately solve the above optimization problem. The algorithm is based on the following steps:

1. Set the initial step size $\alpha = 1$.
2. Compute the residual norm $||r_{k+1}|| = ||F(\textbf{u}_{k+1} - \alpha \Delta \textbf{u})||_2^2$.
3. If the residual norm is smaller than $||r_{k}||$ (the norm of the previous residual), $\alpha$ is accepted and the algorithm terminates. Otherwise, alpha is reduced by a factor `Newton.LineSearchReductionFactor` (default value $0.5$) and we continue with step 2 unless $\alpha$ is smaller or equal the `Newton.LineSearchMinRelaxationFactor` (default value $0.125$) in which we accept the step size, even if the residual norm increased.

@note To customize the line search strategy, create a new class derived from @ref Dumux::NewtonSolver and override the private virtual method @ref Dumux::NewtonSolver::lineSearchUpdate_.
Activate line search update by setting the parameter `Newton.UseLineSearch = true`.

### Chopped update

For instance in problems where the residual has multiple roots, it can be beneficial to restrict the search space to a certain region in the parameter space. This can be achieved by a chopped update strategy, where the solution vector is projected back to a feasible region after each Newton step.
A chopped update strategy can be enabled by setting `Newton.EnableChop = true`. As such strategies are usually problem-dependent, there is no default implementation.  The parameter will merely active the hook @ref Dumux::NewtonSolver::choppedUpdate_ which can be overridden in a derived class, see @ref Dumux::RichardsNewtonSolver for an example implementation using this strategy.

@note To customize the chop update strategy, create a new class derived from @ref Dumux::NewtonSolver and override the private virtual method @ref Dumux::NewtonSolver::choppedUpdate_.
Activate the chopped update strategy by setting the parameter `Newton.EnableChop = true`. The line search and chopped update are mutually exclusive, i.e. if `Newton.UseLineSearch` is set to true, the `Newton.EnableChop` parameter will be ignored.

## Termination criteria

Newton's method finds an approximate solution to the residual equation $F(\textbf{u}_k) \approx \textbf{0}$.
But when is the solution "good enough"? A good criterion makes sure that the residual is small enough but should also be easy to compute.
@ref Dumux::NewtonSolver provides several termination criteria.

The default is the relative shift criterion in addition to bounds on the minimum and maximum number of iterations.
The solver never terminates if $k$ is smaller than `Newton.MinSteps` (default value $2$).
The solver always terminates if $k$ is larger than `Newton.MaxSteps` (default value $18$).
@ref Dumux::NewtonSolver::apply returns a boolean value indicating whether the solver converged or not;
@ref Dumux::NewtonSolver::solve throws @ref Dumux::NumericalProblem is the solver did not converge.

### Relative shift criterion

This criterion is enabled per default. The idea is to terminate the algorithm when the solution vector does not change significantly between iterations.
This criterion does not check the actual residual but is useful in practice and fast to compute. We compute the maximum difference
of any degree of freedom (entry of the vector $\textbf{u}_{k+1}$) between the current and the previous iteration.
This uses the criterion

```math
\text{max}_{i=1,\ldots,n} \left( \frac{|u_{k+1,i} - u_{k,i}|}{\max\lbrace 1.0, |u_{k+1,i} + u_{k,i}| \cdot 0.5 \rbrace} \right) < \epsilon_\mathrm{shift}
```

where $u_{k,i}$ is the $i$-th entry of the vector $\textbf{u}_k$ and
$\epsilon_\mathrm{shift}$ is a user-defined parameter `Newton.MaxRelativeShift` (default value $10^{-8}$).
The relative shift criterion can be disabled by setting `Newton.EnableShiftCriterion = false`.
The relative shift criterion can be combined with a residual criterion (see below) by enabling a residual criterion without disabling the shift criterion.
By default the Newton solver terminates when either of the criteria are fulfilled. If you want to require both criteria to be fulfilled, you can do so
by setting `Newton.SatisfyResidualAndShiftCriterion = true`.

### Relative residual norm criterion

This criterion is disabled by default.
It can be enabled by setting `Newton.EnableResidualCriterion = true`.
The idea is to terminate the algorithm when the norm of the residual vector is small in comparison
to the initial residual norm. The relative residual norm criterion checks if

```math
\frac{||\textbf{r}_{k+1}||_2^2}{||\textbf{r}_{0}||_2^2} < \epsilon_\mathrm{r,rel}^2,
```

where $||\textbf{r}_{k+1}||$ is the norm of the residual vector at iteration $k+1$ and
$||\textbf{r}_{0}||$ is the norm of the residual vector at iteration $0$.
The parameter `Newton.ResidualReduction` (default value $10^{-5}$) defines the threshold $\epsilon_\mathrm{r,rel}$.

### Absolute residual norm criterion

This criterion is disabled by default. The idea is to terminate the algorithm when norm of the residual vector is small enough in comparison with a given absolute threshold.
This requires knowledge about the problem because the magnitude of the residual vector has to be estimated.
The criterion is enabled by setting `Newton.EnableAbsoluteResidualCriterion = true` and  `Newton.EnableResidualCriterion = true`. It checks if

```math
||\textbf{r}_{k+1}||_2^2 < \epsilon_\mathrm{r,abs}^2,
```

where $||\textbf{r}_{k+1}||$ is the norm of the residual vector at iteration $k+1$ and
$\epsilon_\mathrm{r,abs}$ is a user-defined parameter `Newton.MaxAbsoluteResidual` (default value $10^{-5}$).
The absolute residual norm criterion can be combined with a relative shift criterion by enabling both criteria.

## Adaptive time stepping

The extended interface @ref Dumux::NewtonSolver::solve with time loop in the signature provides an automated adaptive time stepping strategy.
The idea is retry the entire algorithm with a lower the time step size if the Newton solver fails to converge. This is repeated until
the Newton solver either converges or the maximum number of allowed retry steps is reached.

The factor by which the time step size is reduced after each non-converged try is set by the parameter `Newton.RetryTimeStepReductionFactor` (default value $0.5$).
The maximum number of retries is set by the parameter `Newton.MaxTimeStepDivisions` (default value $10$). If the maximum number of retries is reached,
the a @ref Dumux::NumericalProblem is thrown.

Finally, @ref Dumux::NewtonSolver::suggestTimeStepSize suggest a new time-step size based on the current time-step size and
the convergence history of the Newton solver. In particular, the implemented heuristic determines
how much the time step size is adjusted depends on how far away the Newton is
from converging within $m$ iterations, where $m$ is set by the parameter `Newton.TargetSteps` (default value $10$).


## Partial reassembly (Quasi-Newton method)

Technically, when using a finite difference approximation of the Jacobian, we use a so-called quasi-Newton method for which
the Jacobian matrix is replaced by an approximation. Another strategy is useful for simulations where the Jacobian matrix is expensive to compute but
not all degrees of freedom are expected to change significantly between iterations (for example if changes only happen in parts of the simulation domain).

If supported by the assembler (e.g. @ref Dumux::FVAssembler), the Newton solver can tell the assembler
to only recompute those entries of the Jacobian matrix that are associated with primary variables that have changed significantly.
This strategy is enabled via the parameter `Newton.EnablePartialReassembly` (default value is `false`).

The threshold for the partial reassembly is computed can also be controlled via parameters.
The threshold $\theta_r$ for each entry is computed as

```math
\theta_r = \max(\theta_\mathrm{min}, \min(\theta_\mathrm{max}, \omega \epsilon_\textrm{shift}))
```
where $\theta_\mathrm{min}$ is set by `Newton.ReassemblyMinThreshold` (default value $10^{-1} \Delta u_\mathrm{shift}$)
and $\theta_\mathrm{max}$ is set by `Newton.ReassemblyMaxThreshold` (default value $10^{-2} \Delta u_\mathrm{shift}$),
$\omega$ is a weight (default value $10^{-3}$), and $\epsilon_\textrm{shift}$
is the maximum relative shift threshold parameter (`Newton.MaxRelativeShift`) described above.

## Primary variable switching

Switching primary variables during Newton's method is a form of nonlinear preconditioning.
We solve an alternative problem with different primary variables where the hope that the root is easier to find.
This can be, for instance, because the finite difference approximation of the Jacobian is numerically better conditioned
in the new variables. The strategy depends on the problem at hand. For an example,
see for instance @ref TwoPNCModel for a model with activated primary variable switch.
