@addtogroup BoxDiscretization

The so called box method is a control volume finite element method (uniting the advantages of the finite-volume (FV) and
finite-element (FE) methods) based on linear nodal basis functions for the trial space.

First, the model domain $\Omega$ is discretized (primary grid, see figure).
Then, a secondary control volume mesh is constructed by connecting the face barycenters and element barycenters,
thus creating a control volume $B_v$, also called box, with node $v$ in the center.
When referring to faces, we mean entities of codimension $1$ with respect to the element dimension.

![Control volume partitioning for the box method.](box.svg){html: width=50%}

Control volumes $B_v$ are partitioned into sub-control volumes (`scv`s) such that each
sub-control volume is the intersection of the control volume with a different primal grid element, $M_v = B_v \cap M$.
The faces of $B_v$ are partitioned into sub-control volume faces $\sigma_k$ (`scvf`s) analogously and $|\sigma_k|$ denotes
the measure of sub-control volume face $k$.
Finally, the integration points $x_k$ which lie on the `scvf`s and the outer normal
vectors $\boldsymbol{n}_{\sigma_k}$ also need to be defined.

The advantage of the FE method is that unstructured grids can be used, while the
FV method is mass conservative. The idea is to apply the FV method (balance of
fluxes across the interfaces) to each control volume and to get the fluxes across
the `scvf`s at the integration points $x_k$ from the FE approach.
Consequently, at each `scvf` the following expression results:

\begin{equation}
   f(\tilde u(x_k)) \cdot \boldsymbol{n}_{\sigma_k} \: |\sigma_k| \qquad \textrm{with}
   \qquad \tilde u(x_k) = \sum_i N_i(x_k) \cdot \hat u_i,
\end{equation}
where $N_i$ represents the basis function of the finite element ansatz at node $i$. The basis functions are defined such that $N_i(x_j)=\delta_{ij}$ with $\delta_{ij}$ being the Kronecker delta.
In the following, the discretization of the balance equation is going to be derived.
From the Reynolds transport theorem follows the general balance equation:

\begin{equation}
  \underbrace{\int_\Omega \frac{\partial}{\partial t} \: u \, \mathrm{d}x}_{1}
  + \underbrace{\int_{\partial\Omega} (\boldsymbol{v} u + \boldsymbol{w}) \cdot \textbf n \, \mathrm{d}\Gamma}_{2} = \underbrace{\int_\Omega q \, \mathrm{d}x}_{3}
\end{equation}

\begin{equation}
  f(u) = \int_\Omega \frac{\partial u}{\partial t} \, \mathrm{d}x + \int_{\Omega} \nabla \cdot
  \underbrace{\left[  \boldsymbol{v} u + \boldsymbol{w}(u)\right] }_{F(u)}  \, \mathrm{d}x - \int_\Omega q \, \mathrm{d}x = 0
\end{equation}
where the first term describes the changes of entity $u$ within a control volume over
time, the second term the advective, diffusive and dispersive fluxes over the interfaces
of the control volume, and the third term is a source or sink. $\Omega$ denotes the
model domain and $F(u) = F(\boldsymbol{v}, p) = F(\boldsymbol{v}(\boldsymbol{x},t), p(\boldsymbol{x},t))$.

Like the finite element method, the box method follows the principle of weighted residuals.
In the function $f(u)$ the unknown $u$ is approximated by discrete values at the
nodes of the primary grid $\hat u_i$ and linear basis functions $N_i$ yielding an
approximate function $f(\tilde u)$.
For instance, an unknown vector $u\in \lbrace \boldsymbol{v}, p, x^\kappa \rbrace$ (velocity, pressure, mole fraction)
is approximated as

\begin{align}
  \tilde p &= \sum_i N_i \hat{p}_i,&
  \tilde{\boldsymbol{v}} &= \sum_i N_i \hat{\boldsymbol{v}}_i,&
  \tilde x^\kappa &= \sum_i N_i \hat x_i^\kappa,&
\end{align}

\begin{align}
  \nabla \tilde p &= \sum_i \nabla N_i \hat{p}_i,&
  \nabla \tilde{\boldsymbol{v}} &= \sum_i \nabla N_i \hat{\boldsymbol{v}}_i,&
  \nabla \tilde x^\kappa  &= \sum_i \nabla N_i \hat x_i^\kappa .&
\end{align}

Due to the approximation in a discrete function space with linear basis functions,
the differential equations are not exactly fulfilled anymore but a residual $\varepsilon$ is produced.

\begin{equation}
  f(u) = 0  \qquad \Rightarrow \qquad f(\tilde u) = \varepsilon
\end{equation}

Application of the principle of weighted residuals, meaning the multiplication
of the residual $\varepsilon$ with a weighting function $W_j$  and claiming that
this product has to vanish within the whole domain,

\begin{equation}
  \int_\Omega \varepsilon W_j \, \mathrm{d}x \overset {!}{=} \: 0 \qquad \textrm{with} \qquad \sum_j W_j =1
\end{equation}
yields the following equation:

\begin{equation}
  \int_\Omega \frac{\partial \tilde u}{\partial t} W_j \, \mathrm{d}x + \int_\Omega
   \left[ \nabla \cdot F(\tilde u) \right] W_j  \, \mathrm{d}x - \int_\Omega q W_j \, \mathrm{d}x = \int_\Omega \varepsilon W_j \, \mathrm{d}x \: \overset {!}{=} \: 0.
\end{equation}

For Galerkin schemes (standard finite element method), the weighting functions $W_j$ are chosen the same as the ansatz functions $N_j$.
However, this does not yield a locally mass-conservative discretization scheme.
Instead, for the Box method, the weighting functions $W_j$ are chosen as the piece-wise constant functions over a
control volume box $B_j$,
\begin{equation}
  W_j(x) = \begin{cases}
            1 &x \in B_j \\
      0 &x \notin B_j.\\
           \end{cases}
\end{equation}
Thus, the Box method is a Petrov-Galerkin scheme with weighting functions
not belonging to the same function space as the ansatz functions.

Inserting the introduced weighting functions and using the divergence theorem results in
\begin{equation}
  \int_{B_j} \frac{\partial \tilde u}{\partial t} \, \mathrm{d}x + \int_{\partial B_j}  F(\tilde u) \cdot \boldsymbol{n} \, \mathrm{d}\Gamma_{B_j} - \int_{B_j} q \, \mathrm{d}x  \overset {!}{=} \: 0,
\end{equation}
which has to hold for every control volume $B_j$.

The first term in previous equation can be written as
\begin{equation}
\int_{B_j} \frac{\partial \tilde u}{\partial t} \, \mathrm{d}x = \frac{d}{dt} \int_{B_j} \sum_i \hat u_i N_i  \, \mathrm{d}x = \sum_i \frac{\partial \hat u_i}{\partial t} \int_{B_j}  N_i  \, \mathrm{d}x.
\end{equation}
A technique called mass lumping is applied by assuming that the storage capacity is
reduced to the nodes. This means that the integrals $M_{i,j} = \int_{B_j}  N_i \, \mathrm{d}x$
are replaced by some mass lumped terms $M^{lump}_{i,j}$ which are defined as
\begin{equation}
   M^{lump}_{i,j} =
   \begin{cases}
        |B_j| &j = i\\
        0 &j \neq i,\\
    \end{cases}
\end{equation}
where $|B_j|$ is the volume of the control volume $B_j$ associated with node $j$.
Mass lumping yields
\begin{equation}
  |B_j| \frac{\partial \hat u_j}{\partial t}
  +  \int_{\partial B_j}  F(\tilde u) \cdot \boldsymbol{n} \, \mathrm{d}\Gamma_{B_j} - Q_j = 0,
\end{equation}
where $Q_j$ is an approximation (using some quadrature rule) of the integrated source/sink term $\int_{B_j} q \, \mathrm{d}x$.

Using an implicit Euler time discretization finally
leads to the discrete form
\begin{equation}
  |B_j| \frac{\hat u_j^{n+1} - \hat u_j^{n}}{\Delta t}
  + \int_{\partial B_j}  F(\tilde u^{n+1}) \cdot \boldsymbol{n}
  \;   \mathrm{d}\Gamma_{B_j} - Q_j^{n+1} \: = 0,
\end{equation}
which is to be fulfilled for each box $B_j$.
