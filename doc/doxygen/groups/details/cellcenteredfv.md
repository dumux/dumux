@addtogroup CCDiscretization

Cell-centered finite volume methods use the elements of the grid as control volumes.
For each control volume the discrete values are determined at the element/control
volume center (not required to be the barycenters but often chosen so).

The main ideas of the two-point flux approximation and the multi-point flux approximation methods are described in
the documentation pages for @ref CCTpfaDiscretization and @ref CCMpfaDiscretization.

Please also note that other types of equations, e.g. instationary parabolic problems, can be discretized by applying some time discretization scheme to the time derivatives and by using the finite-volume scheme for the flux discretization. For simplicity the discussion here is restricted to the elliptic problem presented above.
