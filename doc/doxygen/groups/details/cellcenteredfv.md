@addtogroup CCDiscretization

Cell-centered finite volume methods use the elements of the grid as control volumes.
For each control volume the discrete values are determined at the element/control
volume center (not required to be the barycenters but often chosen so).

We consider a domain $\Omega \subset \mathbb{R}^d$, $d \in \{ 2, 3 \}$ with boundary $\Gamma = \partial \Omega = \Gamma_D \cup \Gamma_N$, and the following elliptic problem with Dirichlet boundary conditions on $\Gamma_D$ and Neumann boundary on $\Gamma_N$:
\begin{equation}
  \begin{aligned}
                   \nabla \cdot \left( - \boldsymbol{\Lambda} \nabla u \right) &= q   &&\mathrm{in} \, \Omega \\
               \left( - \boldsymbol{\Lambda} \nabla u \right) \cdot \boldsymbol{n} &= v_N &&\mathrm{on} \, \Gamma_N \\
                                                                   u &= u_D &&\mathrm{on} \, \Gamma_D.
    \label{eq:elliptic}
  \end{aligned}
\end{equation}

Here, $\boldsymbol{\Lambda} = \boldsymbol{\Lambda}(\boldsymbol{x}, \boldsymbol{u})$ is a symmetric and positive definite tensor of second rank (e.g. permeability, diffusivity, etc.), $u = u (\boldsymbol{x})$ is unknown and $q = q(\boldsymbol{x}, \boldsymbol{u})$ is a source/sink.
We denote by $\mathcal{M}$ the mesh that results from the division of the domain $\Omega$ into $n_e$ control volumes $K \subset \Omega$. Each $K$ is a polygonal open set such that $K \cap L = \emptyset, \forall{K \neq L}$ and $\overline{\Omega} = \cup_{K \in \mathcal{M}} \overline{K}$.

For the derivation of the finite volume formulation, we integrate the governing equation over a control volume $K$ and apply the Gauss divergence theorem,

\begin{equation}
    \int_{K} \nabla \cdot \left( - \boldsymbol{\Lambda} \nabla u \right) \, \mathrm{d} \Omega = \int_{\partial K} \left( - \boldsymbol{\Lambda} \nabla u \right) \cdot \boldsymbol{n} \, \mathrm{d} \Gamma = \int_K q \, \mathrm{d}x.
    \label{eq:ellipticIntegrated}
\end{equation}

Splitting the control volume boundary $\partial K$ into a finite number of faces $\sigma \subset \partial K$ (such that $\sigma = \overline{K} \cap \overline{L}$ for some neighboring control volume $L$) and replacing the exact fluxes by an approximation, i.e. $F_{K, \sigma} \approx \int_{\sigma} \left( - \boldsymbol{\Lambda}_K \nabla u \right) \cdot \boldsymbol{n} \mathrm{d} \Gamma$ (here $\boldsymbol{\Lambda}_K$ is the value of $\boldsymbol{\Lambda}$ associated with control volume $K$), yield the discrete equation
\begin{equation}
    \sum_{\sigma \subset \partial K} F_{K, \sigma} = Q_K, \quad \forall \, {K \in \mathcal{M}},
\label{eq:ccdisc}
\end{equation}
where $F_{K, \sigma}$ is the discrete flux through face $\sigma$ flowing out of cell $K$ and $Q_K := \int_K q \, \mathrm{d}x$ is the integrated source/sink term. This is the prototypical cell-centered finite volume formulation.
Finite volume schemes differ in the way how the term
$(\boldsymbol{\Lambda}_K \nabla u ) \cdot \boldsymbol{n} $ is approximated (i.e. the choice of the fluxes $F_{K, \sigma}$). Using the symmetry of the tensor $\boldsymbol{\Lambda}_K$, the flux can be rewritten as
$\nabla u  \cdot \boldsymbol{\Lambda}_K\boldsymbol{n}$, which corresponds to the directional derivative of $u$ in co-normal direction $\boldsymbol{\Lambda}_K\boldsymbol{n}$.

The main ideas of the two-point flux approximation and the multi-point flux approximation methods are described in
the documentation pages for @ref CCTpfaDiscretization and @ref CCMpfaDiscretization.

Please also note that other types of equations, e.g. instationary parabolic problems, can be discretized by applying some time discretization scheme to the time derivatives and by using the finite-volume scheme for the flux discretization. For simplicity the discussion here is restricted to the elliptic problem presented above.
