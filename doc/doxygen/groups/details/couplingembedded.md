@addtogroup EmbeddedCoupling

![Examples for possible systems with embedded coupling.](multidomainEmbeddedCoupling.svg){html: width=55%}

Typical examples include embedded one-dimensional networks for the simulation of blood tissue perfusion or root-soil interaction, and also embedded fracture models. The latter ones usually embed two-dimensional planes in a three-dimensional bulk and allow simulating fractures in porous media. Most of the time, in DuMu<sup>x</sup> there is no extra coupling mapper that adopts the management of coupling stencils. Instead, the coupling manager incorporates this subtask. The general procedure for embedded coupling starts by intersecting the given grids of differing dimensions. Then, quadrature points are placed along the intersections which serve as point sources. Finally, for each intersection a source term can be computed by integrating over all quadrature points which hold individual integration weights and integration elements.
