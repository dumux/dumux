# Basic numerics

In this section we discuss some basic numerical concepts for the
solution of partial differential equations used in DuMux.

## Newton's method

Nonlinear systems of equations can be solved with Newton's method.
Newton's method is quadratically convergent, if the initial solution is
close enough to the root of the residual equation. It does however not
guarantee global convergence.

The discrete differential equations are formulated in residual form. All terms are
on the left hand side and are summed up. The terms contain values for the primary
variables which are part of the solution vector $\textbf{u}$. The sum of the terms
is called residual $\textbf{r}(\textbf{u})$ which is a function of the solution. For
example:

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

We are looking for the unknown solution $\textbf{u}$ and use Newton's method to
create a series of $\textbf{u}^i$ that we hope converges to $\textbf{u}$.
We start with an initial guess $\textbf{u}^0$ and
calculate the initial residual $\textbf{r}(\textbf{u}^0)$. Next,
we calculate the derivative of the residual with respect to the solution.
This is the Jacobian matrix

\begin{align*}
  \frac{\text{d}}{\text{d}\textbf{u}}\textbf{r} \left(\textbf{u}^i\right)
  = J_{\textbf{r} \left(\textbf{u}^i\right)}
  = \left(\frac{\text{d}}{\text{d}\textbf{u}^i_m}\textbf{r} \left(\textbf{u}^i\right)_n\right)_{m,n}
\end{align*}

with $i$ denoting the Newton iteration step.
Each column is the residual derived with respect to the $m$th entry of $\textbf{u}^i$.

The Jacobian indicates the direction where the residual increases. By solving the
linear system

\begin{align*}
  J_{\textbf{r}(\textbf{u}^i)} \cdot \textbf{x}^i = \textbf{r}(\textbf{u}^i)
\end{align*}

we calculate the direction of maximum growth $\textbf{x}^i$. We subtract it from
our current solution to get a new, better solution
$\textbf{u}^{i+1} = \textbf{u}^i - \textbf{x}^i$.

We repeat the calculation of the Jacobian $J_{\textbf{r}(\textbf{u}^i)}$ and the
direction of maximum growth $\textbf{x}^i$ until our approximated solution becomes good enough.
See `Dumux::NewtonSolver` for the implementation of the Newton method based solver in DuMux
which features various convergence criteria and a simple line search algorithm to improve the update.

## Time discretization

In the following, we describe two basic time stepping schemes which are first-order accurate.
Our systems of partial differential equations are discretized in space and in time.

Let us consider the general case of a balance equation of the following form

\begin{equation}
\frac{\partial m(u)}{\partial t} + \nabla\cdot\boldsymbol{f}(u, \nabla u) + q(u) = 0,
\end{equation}

seeking an unknown quantity $u$ in terms of storage $m$, flux $\boldsymbol{f}$ and source $q$.
All available Dumux models can be written mathematically in this form
with possibly vector-valued quantities $u$, $m$, $q$ and a tensor-valued flux $\boldsymbol{f}$.
For the sake of simplicity, we here assume scalar quantities $u$, $m$, $q$ and a vector-valued
flux $\boldsymbol{f}$.

To discretize the continuous form of the balance equations, we need to choose an
approximation for the temporal derivative $\partial m(u)/\partial t$.
While many elaborate methods for this approximation exist,
we focus on the simplest one of a first order difference quotient

\begin{equation}
\frac{\partial m(u_{k/k+1})}{\partial t}
\approx \frac{m(u_{k+1}) - m(u_k)}{\Delta t_{k+1}}
\end{equation}

for approximating the solution $u$ at time $t_k$ (forward) or $t_{k+1}$ (backward).
The question of whether to choose the forward or the backward quotient leads to the
explicit and implicit Euler method, respectively.
Using the explicit Euler method leads to

\begin{equation}
\frac{m(u_{k+1}) - m(u_k)}{\Delta t_{k+1}} + \nabla\cdot\boldsymbol{f}(u_k, \nabla u_k) + q(u_k) = 0,
\end{equation}

whereas the implicit Euler method is described as

\begin{equation}
\frac{m(u_{k+1}) - m(u_k)}{\Delta t_{k+1}}
+ \nabla\cdot\boldsymbol{f}(u_{k+1}, \nabla u_{k+1}) + q(u_{k+1}) = 0.
\end{equation}

Once the solution $u_k$ at time $t_k$ is known, it is straightforward
to determine $m(u_{k+1})$,
while attempting to do the same for the equation discreized with the implicit Euler
involves the solution of a system of equations.
The explicit method is stable only
if the time step size $\Delta t_{k+1}$ is below a certain limit that depends
on the specific balance equation, whereas the implicit method
is unconditionally stable.
