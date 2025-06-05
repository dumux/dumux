\addtogroup FacetCoupling

In this approach, a bulk subproblem of dimension d and a lower-dimensional subproblem of dimension d-1 are coupled. The lower-dimensional subproblem is embedded on the facets of the bulk. For facet coupling, an appropriate mesh file that contains the bulk elements but also the lower-dimensional elements needs to be provided. DuMu<sup>x</sup> then extracts a coupling map that contains neighborhood information on the coupling elements in both grids.

![Example for conforming mixed-dimension facet coupling.](multidomainFacetCoupling.svg){html: width=55%}

Possible examples span discrete facet conforming fracture models and problems with physics on a domain surface.
