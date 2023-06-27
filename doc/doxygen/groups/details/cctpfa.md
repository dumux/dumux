@addtogroup CCTpfaDiscretization

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

![Two neighboring control volumes, K, L, sharing the face σ](cctpfa.svg){html: width=50%}

The linear two-point flux approximation is a simple but robust cell-centered finite-volume scheme, which is commonly used in commercial software. The scheme can be derived by using the co-normal decomposition, which reads
\begin{equation}
\boldsymbol{\Lambda}_K \boldsymbol{n}_{K, \sigma} = t_{K,\sigma} \boldsymbol{d}_{K,\sigma} + \boldsymbol{d}^{\bot}_{K,\sigma}, \quad  t_{K,\sigma} = \frac{\boldsymbol{n}_{K, \sigma}^T \boldsymbol{\Lambda}_K \boldsymbol{d}_{K,\sigma} }{\boldsymbol{d}_{K,\sigma}^T \boldsymbol{d}_{K,\sigma}}, \; \boldsymbol{d}^{\bot}_{K,\sigma} = \boldsymbol{\Lambda}_K \boldsymbol{n}_{K, \sigma} - t_{K,\sigma} \boldsymbol{d}_{K,\sigma},
\label{eq:conormalDecTpfa}
\end{equation}
with the tensor $\boldsymbol{\Lambda}_K$ associated with control volume $K$, the distance vector $\boldsymbol{d}_{K,\sigma} := \boldsymbol{x}_\sigma - \boldsymbol{x}_K$ and $\boldsymbol{d}_{K,\sigma}^T \boldsymbol{d}^{\bot}_{K,\sigma} = 0$ (orthogonality). The same can be done for the conormal $\boldsymbol{\Lambda}_L \boldsymbol{n}_{L, \sigma}$. The $t_{K,\sigma}$ and $t_{L,\sigma}$ are the transmissibilities associated with the face $\sigma$. These transmissibilities are calculated in DuMux with the free function Dumux::computeTpfaTransmissibility.

With the introduced notation, it follows that for each cell $K$ and face $\sigma$
\begin{equation}
\nabla u \cdot \boldsymbol{\Lambda}_K \boldsymbol{n}_{K, \sigma} =  t_{K,\sigma} \nabla u \cdot \boldsymbol{d}_{K,\sigma} + \nabla u \cdot \boldsymbol{d}^{\bot}_{K,\sigma}.
\end{equation}
For the Tpfa scheme, the second part in the above equation is neglected. By using the fact that $\nabla u \cdot \boldsymbol{d}_{K,\sigma} \approx u_\sigma - u_K$, the discrete fluxes for face $\sigma$ are given by
\begin{equation}
F_{K,\sigma} = -\vert{\sigma}\vert  t_{K,\sigma} (u_\sigma - u_K), \qquad F_{L,\sigma} = -\vert{\sigma}\vert  t_{L,\sigma} (u_\sigma - u_L).
\end{equation}
Enforcing local flux conservation, i.e. $F_{K,\sigma}+F_{L,\sigma}=0$, results in
\begin{equation}
u_\sigma = \frac{t_{K,\sigma} u_K + t_{L,\sigma} u_L}{t_{K,\sigma}  + t_{L,\sigma}}.
\end{equation}
With this, the fluxes $F_{K,\sigma}$ are rewritten as
\begin{equation}
F_{K,\sigma} = \vert{\sigma}\vert \frac{t_{K,\sigma} t_{L,\sigma}}{t_{K,\sigma} + t_{L,\sigma}} (u_K - u_L), \quad F_{L,\sigma} = \vert{\sigma}\vert  \frac{t_{K,\sigma} t_{L,\sigma}}{t_{K,\sigma} + t_{L,\sigma}} (u_L - u_K).
\label{eq:TPFAFlux}
\end{equation}
By neglecting the orthogonal term, the consistency of the scheme is lost for general grids, where $\nabla u \cdot \boldsymbol{d}^{\bot}_{K,\sigma} \not = 0$. The consistency is achieved only for so-called K-orthogonal grids for which $\boldsymbol{d}^{\bot}_{K,\sigma} = 0$. For such grids we deduce that
\begin{equation}
\frac{t_{K,\sigma} t_{L,\sigma}}{t_{K,\sigma} + t_{L,\sigma}} = \frac{\tau_{K,\sigma} \tau_{L,\sigma}}{\tau_{K,\sigma} d_{L,\sigma} + \tau_{L,\sigma} d_{K,\sigma}},
\end{equation}
with $\tau_{K,\sigma} := \boldsymbol{n}_{K, \sigma} \boldsymbol{\Lambda}_K\boldsymbol{n}_{K, \sigma}, \tau_{L,\sigma} := \boldsymbol{n}_{L, \sigma} \boldsymbol{\Lambda}_L\boldsymbol{n}_{L, \sigma}$, $d_{K,\sigma}:= \boldsymbol{n}_{K, \sigma} \cdot \boldsymbol{d}_{K, \sigma}$, and $d_{L,\sigma}:= \boldsymbol{n}_{L, \sigma} \cdot \boldsymbol{d}_{L, \sigma}$. This reduces, for the case of scalar permeability, to a distance weighted harmonic averaging of permeabilities.
