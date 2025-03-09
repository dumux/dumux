@addtogroup Nonlinear

The Newton solver is used to solve a non-linear system of equations.

The discrete differential equations are formulated in residual form.
This means that all terms of each equation are brought to one side and summed up, which forms a component of the residual vector $\textbf{r}(\textbf{u})$, which is a function of the solution vector $\textbf{u}$. The goal is to find the solution vector $\textbf{u}$ that makes the residual vector equal to zero (or very close to zero).

For example:

\begin{align*}
\underbrace{
  \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left(
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \boldsymbol{K}
 \left(\nabla p_\alpha - \varrho_{\alpha} \boldsymbol{g} \right)
 \right) - q_\alpha} _
{=: \, \textbf{r}(\textbf{u})}
= 0
\end{align*}

Since Newton's method is an iterative algorithm, meaning it generates a sequence of approximations, $\textbf{u}^i$, to converge to the solution vector, that satisfies $\textbf{r}(\textbf{u}) = 0$.
We start with an initial guess $\textbf{u}^0$ and
calculate the initial residual $\textbf{r}(\textbf{u}^0)$. Then,
we calculate the derivative of the residual with respect to the solution, which gives the Jacobian matrix, $J_{\textbf{r} \left(\textbf{u}^i\right)}$, representing the sensitivity of the residual to changes in the solution. It linearizes the nonlinear system around the current approximation $\textbf{u}^i$:

\begin{align*}
  \frac{\text{d}}{\text{d}\textbf{u}}\textbf{r} \left(\textbf{u}^i\right)
  = J_{\textbf{r} \left(\textbf{u}^i\right)}
  = \left(\frac{\text{d}}{\text{d}\textbf{u}^i_m}\textbf{r} \left(\textbf{u}^i\right)_n\right)_{m,n}
\end{align*}

By solving the linear system

\begin{align*}
  J_{\textbf{r}(\textbf{u}^i)} \cdot \textbf{x}^i = \textbf{r}(\textbf{u}^i)
\end{align*}

we calculate the direction of maximum growth $\textbf{x}^i$ and subtract it from
our current approximation to get a new, better approximation
$\textbf{u}^{i+1} = \textbf{u}^i - \textbf{x}^i$.

We repeat the calculation of the Jacobian $J_{\textbf{r}(\textbf{u}^i)}$ and the
direction of maximum growth $\textbf{x}^i$ until our approximated solution becomes good enough.

See `Dumux::NewtonSolver` for the implementation of the Newton method based solver in DuMux which features various convergence criteria and a simple line search algorithm to improve the update.


 Users have the flexibility to adjust various Newton solver parameters to optimize convergence and efficiency.
 The parameters, that can be set in the input file, e.g., `params.input`, are:
 - **UseLineSearch**: the default value is false. If is set to true, an inexact line search method is used, which multiplies the full Newton step predicted by Newton method with a relaxation factor, $\alpha$ to improve the convergence behavior, i.e., increase the residual reduction.
  \begin{align*}
  \textbf{u}^{i+1} = \textbf{u}^i - \alpha\textbf{x}^i
  \end{align*}
  The relaxation factor decreases repeatedly until the residual reduction does not improve anymore or the relaxation factor reaches a minimum value, i.e., `LineSearchMinRelaxationFactor`.
 - **LineSearchMinRelaxationFactor**: The default value is 0.125. to set the minimum values of the line search relaxation factor.
 - **EnableChop**: the default value is false. If is set to true, the Newton step is reduced. Please note that both `UseLineSearch` and `EnableChop` can be set to true simultaneously. To this end, the chopped Newton solver strategy is not yet implemented. Therefore, even if `EnableChop` is set to true, the Newton step will not be reduced in the current version.
 - **EnableShiftCriterion**: the default value is true. The Newton solver converges, when the maximum relative difference between the solution vector from the last iteration and the one before that is less than the `MaxRelativeShift`.
 - **EnableResidualCriterion**: the default value is false. The Newton converges, when the relative reduction of the residual vector to the initial residual vector, i.e., at the beginning of the time-step, is less than `ResidualReduction` (default residual criterion). When `EnableAbsoluteResidualCriterion` is also set to true, then the Newton converges, when the norm of the residual vector is less than `MaxAbsoluteResidual`. It should be noted that setting the `EnableResidualCriterion = true` **does not** automatically set `EnableShiftCriterion` to false. That means, if only `EnableResidualCriterion` is set to true, then the Newton converges when either the residual criterion or the shift criterion is fulfilled (the shift criterion is checked first). To use only the residual criterion, 'EnableShiftCriterion' must explicitly be set to false.
 - **EnableAbsoluteResidualCriterion**: the default value is false. By setting it to true, the residual criterion changes to when the norm of the residual vector is less than `MaxAbsoluteResidual`.
 - **SatisfyResidualAndShiftCriterion**: the default value is false. The Newton solver converges when both residual criterion and shift criterion are fulfilled.
 - **MaxRelativeShift**: the default value is 1e-8.
 - **MaxAbsoluteResidual**: the default value is 1e-5.
 - **ResidualReduction**: the default value is 1e-5.
 - **MinSteps**: the default value is 2. The minimum iterations, that Newton must do before it converges.
 - **MaxSteps**: the default value is 18. The maximum iterations, that Newton is allowed to do until it converges. If Newton can not converges after doing this number of iterations, the last maximum relative shift in the solution vector or the reduction in the residual is checked. If the last maximum relative shift was reduced by a factor of at least 4 in comparison to the one before it, Newton is allowed to proceed to the next iteration. Otherwise the time-step will be redone with a reduced time-step size.
 - **TargetSteps**: the default value is 10. It sets the optimum number of iterations we want to achieve. It plays a role in time-step size suggestion by Newton for the next time-step, considering how fast Newton converged, i.e., the number of iterations done to converge in a time-step (`numSteps_`). That means, if the current `numSteps_` is below `TargetSteps`, Newton will suggest a larger time-step size for the next time step, and a smaller time-step size if it's above.
 - **RetryTimeStepReductionFactor**: the default value is 0.5. If Newton does not converges, then the time-step size will be reduced by a factor of `RetryTimeStepReductionFactor` and the time-step will be redone with the new time-step size.
 - **MaxTimeStepDivisions**: the default value is 10. If even after `MaxTimeStepDivisions` of time-step size reduction, the Newton solver does not succeed to solve the system, the simulation will be aborted.
 - **EnablePartialReassembly**: the default value is false. If it is set to true, the Jacobian matrix is partially reassembled in the next iteration, computing only components associated with primary variables that have changed beyond the Reassembly threshold. The threshold is computed as: <br> Reassembly threshold = max(`ReassemblyMinThreshold`, min(`ReassemblyMaxThreshold`, `ReassemblyShiftWeight` $\times$ Relative shift))

 - **ReassemblyMinThreshold**: the default value is $10^{-1} \times$ `MaxRelativeShift`.
 - **ReassemblyMaxThreshold**: the default value is $10^2 \times$ `MaxRelativeShift`.
 - **ReassemblyShiftWeight**: the default value is $10^{-3}$.
