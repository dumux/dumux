@addtogroup DiamondDiscretization

The diamond method is a control-volume finite element scheme with piecewise linear discontinuous basis functions
(Crouzeix-Raviart on simplices and Rannacher-Turek on quadrilaterals). The basis functions are continuous only at
the element face centers but are generally discontinuous along the facet and thus also at vertices.
Control volumes are constructed around the face centers of the element faces of the primary grid.
This results in diamond-shaped control volumes that form a dual grid.

![Control volume partitioning for the diamond CVFE method.](fcdiamond.svg){html: width=50%}

The diamond CVFE scheme can, for example, be used for the velocity combined with cell-centered FV schemes for the pressure
to discretize the incompressible Stokes equation (in formulation with full velocity gradient in the flux term; for the formulation
with symmetric velocity gradient additional face stabilization terms are required) on unstructured grids.
