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

Since Newton's method is an iterative algorithm, meaning it generates a sequence of approximations, $\textbf{u}^i$, to converge to the solution vector that satisfy $\textbf{r}(\textbf{u}) = 0$.
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
 The aparmeters, that can be set in the input file, e.g., `params.input`, are:
 * UseLineSearch: if is set to true, an inexact line search method is used, which multiplies the full Newton step predicted by Newton method with a relaxation factor to improve the convergence behavior, i.e., increase the residual reduction. The relaxation factor decreases repeatedly until the residual reduction does not increase anymore or it reaches a minimum value.
 * LineSearchMinRelaxationFactor: to set the minimum values of the line search relaxation factor. The default value is 0.125.
 * EnableChop:
 * EnablePartialReassembly:
 * EnableAbsoluteResidualCriterion:
 * EnableShiftCriterion:
 * EnableResidualCriterion:
 * SatisfyResidualAndShiftCriterion:
 * MaxRelativeShift:
 * MaxAbsoluteResidual:
 * ResidualReduction:
 * MinSteps:
 * MaxSteps:
 * TargetSteps:
 * ReassemblyMinThreshold:
 * ReassemblyMaxThreshold:
 * ReassemblyShiftWeight:
 * RetryTimeStepReductionFactor:
 * MaxTimeStepDivisions:
