@addtogroup CCMpfaDiscretization

We continue from the introduction to @ref CCDiscretization using as an example the following scalar elliptic
problem for the unknown $u$,
\begin{equation}
  \begin{aligned}
                   \nabla \cdot \left( - \boldsymbol{\Lambda} \nabla u \right) &= q   &&\mathrm{in} \, \Omega \\
               \left( - \boldsymbol{\Lambda} \nabla u \right) \cdot \boldsymbol{n} &= v_N &&\mathrm{on} \, \Gamma_N \\
                                                                   u &= u_D &&\mathrm{on} \, \Gamma_D,
    \label{eq:elliptic}
  \end{aligned}
\end{equation}
and recall that after integration of the governing equation over a control volume $K$ and application of
the divergence theorem, we replace exact fluxes by numerical approximations $F_{K, \sigma} \approx \int_{\sigma} \left( - \boldsymbol{\Lambda}_K \nabla u \right) \cdot \boldsymbol{n} \mathrm{d} \Gamma$.

For the Mpfa method, expressions for the face fluxes $F_{K, \sigma}$ are obtained by introducing intermediate face unknowns $u_\sigma$ in addition to the cell unknowns $u_K$ and enforcing the physically motivated continuity of fluxes and continuity of the solution across the faces. For a face $\sigma$ between the two polygons $K$ and $L$ these conditions read:
\begin{equation}
    \begin{aligned}
        &F_{K, \sigma} + F_{L, \sigma} = 0 \\
        &{u}_{K,\sigma} = {u}_{L,\sigma} = {u}_{\sigma}.
        \label{eq:sigmaConditions}
    \end{aligned}
\end{equation}
Using these conditions, the intermediate face unknowns ${u}_\sigma$ can be eliminated and the fluxes are expressed as a function of the cell unknowns $u_N$ and associated transmissibilities $t^N_{K,\sigma}$:

\begin{equation}
    F_{K,\sigma} = \sum_{N \in \mathcal{S}_{K,\sigma}} t^N_{K,\sigma} u_{N}.
    \label{eq:FVFluxExpression}
\end{equation}

![Interaction region for the Mpfa-O method. Three elements belong to the interaction region around vertex v. The fluxes over the sub-control volume faces therefore depend on the three cell unknowns uK, uL, uM.](ccmpfa.svg){html: width=50%}

The main difference between different finite volume schemes is the choice of the flux expression, i.e. the computation of the $t^N_{K,\sigma}$ and the size of $\mathcal{S}_{K,\sigma}$. For the Tpfa method this is presented in @ref CCTpfaDiscretization. There, the stencil and transmissibilities are given as
\begin{equation*}
\mathcal{S}_{K,\sigma} = \lbrace K,L \rbrace, \quad t^K_{K,\sigma} =  \vert{\sigma}\vert  \frac{t_{K,\sigma} t_{L,\sigma}}{t_{K,\sigma} + t_{L,\sigma}},\; t^L_{K,\sigma} =  -\vert{\sigma}\vert  \frac{t_{K,\sigma} t_{L,\sigma}}{t_{K,\sigma} + t_{L,\sigma}},
\end{equation*}
with the transmissibilities $t_{K,\sigma}, t_{L,\sigma}$ (see @ref CCTpfaDiscretization).

In the following, we present a multi-point flux approximation method called Mpfa-O method introduced in @cite A3:aavatsmark:2002. The main difference to the Tpfa scheme is the fact that a consistent discrete gradient is constructed, i.e. the term $\nabla u \cdot \mathbf{d}^{\bot}_{K,\sigma}$ is not neglected.

For this scheme, a dual grid is created by connecting the barycenters of the cells with the barycenters of the faces ($d=2$) or the barycenters of the faces and edges ($d=3$).
This divides each primary grid face into $n$ sub-control volume faces $\sigma^v$, where $n$ is the number of corners of the primary grid face and the superscript $v$ refers to the vertex the
sub-face can be associated with (see figure above). Also, it sub-divides the control volume $K$ into regions $K_v$ that can be associated with the vertex $v$. These regions are used to construct a local interpolation scheme.
We allow for piecewise constant $\mathbf{\Lambda}$ (denoted as $\mathbf{\Lambda}_K$ for each cell $K$)
and construct discrete gradients $\nabla_\mathcal{D}^{K_v} u$ in each region $K_v$.
In the following, we restrict our discussion to the two-dimensional setup that is shown in the figure above.
Here, the discrete gradients are constructed to be consistent such that the following conditions hold:
\begin{equation}
\nabla_\mathcal{D}^{K_v} u \cdot (\mathbf{x}_{\sigma^v_1}- \mathbf{x}_{K}) = u_{\sigma^v_1} - u_K, \quad \nabla_\mathcal{D}^{K_v} u \cdot (\mathbf{x}_{\sigma^v_3}- \mathbf{x}_{K}) = u_{\sigma^v_3} - u_K.
\end{equation}
Thus, a discrete gradient in $K_v$ that fulfills these conditions is given as
\begin{equation}
\nabla_\mathcal{D}^{K_v} u  = \mathbb{D}^{-T}_{K_v}
 \begin{bmatrix}
  u_{\sigma^v_1} - u_K \\
  u_{\sigma^v_3} - u_K
 \end{bmatrix}, \qquad \text{ with }\; \mathbb{D}_{K_v} :=
  \begin{bmatrix}
   \mathbf{x}_{\sigma^v_1}- \mathbf{x}_K & \mathbf{x}_{\sigma^v_3} - \mathbf{x}_K
 \end{bmatrix}.
\end{equation}

This enables us to write the discrete flux across $\sigma^v_1$ from cell $K$ as follows:
\begin{equation}
    F_{K, \sigma^v_1} := - |\sigma^v_1| \mathbf{n}_{\sigma^v_1}^T \mathbf{\Lambda}_K \nabla_\mathcal{D}^{K_v} u.
    \label{eq:discreteFlux}
\end{equation}
Inserting the discrete gradient yields
\begin{equation}
    F_{K, \sigma^v_1} = \omega_{K,\sigma^v_1\sigma^v_1}(u_K - u_{\sigma^v_1}) + \omega_{K,\sigma^v_1 \sigma^v_3}(u_K - u_{\sigma^v_3}),
    \label{eq:discreteFluxRef}
\end{equation}
with $(\omega_{K,\sigma^v_1\sigma^v_1},\omega_{K,\sigma^v_1 \sigma^v_3})^T = |\sigma^v_1| \mathbb{D}^{-1}_{K_v}\mathbf{\Lambda}_K \mathbf{n}_{\sigma^v_1}$. These values are calculated in DuMux by using the function Dumux::computeMpfaTransmissibility.


To deduce a cell-centered scheme, the introduced face unknowns $u_{\sigma^v_i}$ have to be eliminated. This is done by enforcing flux continuity for each sub-control volume face, i.e.
\begin{align}
F_{K, \sigma^v_1} + F_{L, \sigma^v_1} &= 0, \\ F_{K, \sigma^v_3} + F_{M, \sigma^v_3} &= 0, \\ F_{L, \sigma^v_2} + F_{M, \sigma^v_2} &= 0.
\end{align}
This results in a system of equations for the face unknowns $\mathbf{u}_{\sigma}$
\begin{equation}
\mathbb{A}^{3\times 3} \mathbf{u}_{\sigma} = \mathbb{B}^{3\times 3} \mathbf{u},
\end{equation}
where $\mathbf{u}$ contains the three cell unknowns $u_K,u_L,u_M$ and $\mathbf{u}_{\sigma}$ the three face unknowns $u_{\sigma^v_1}, u_{\sigma^v_2}, u_{\sigma^v_3}$.
Inserting these face unknowns into the flux expression yields
\begin{equation}
    F_{K,\sigma^v_i} = \sum_{N \in \lbrace K,L,M \rbrace } t^N_{K,\sigma^v_i} u_{N} = \mathbf{t}_{K,\sigma^v_i} \cdot \mathbf{u},
\end{equation}
for each cell $K$ and sub-control volume face $\sigma^v_i$.
