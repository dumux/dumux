@addtogroup FaceCenteredStaggeredDiscretization

The finite volume method on a staggered grid is a classical discretization scheme for the momentum balance of the Navier-Stokes equations  \cite Harlow1965.
It is constructed on Cartesian (tensor product) grids. Here, we use the same abstraction as for other finite volume and CVFE schemes and split
the staggered (shifted) control volumes into element-wise parts: sub control volumes. Fluxes are assembled over sub control volumes faces,
where two control volume faces are fully contained in an element (these are called frontal faces), and two faces are shared between
neighboring elements and therefore split into two sub control volume faces (these are called lateral faces).

![Control volume partitioning for the staggered method.](fcstaggered.svg){html: width=50%}

The image shows the control volume partition in 2D.
It also shows the neighboring degrees of freedom involved in the interpolation of gradients at the flux integration points.
