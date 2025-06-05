@addtogroup BoundaryCoupling

The boundary coupling mode is excellent for coupling subproblems that touch at domain boundaries which also means that each subproblem owns an individual domain. The subproblems may have identical or different dimensionalities, but can also consist of identical or differing physical models. Due to potential differences in DOF location of the the different models, the coupling manager needs to assure that the subproblems' data is mapped correctly to each other at the coupling boundary. This is particularly important when the subproblems use different discretizations and the degrees of freedom do not coincide at the same positions.

![Examples for possible systems with boundary coupling.](multidomainBoundaryCoupling.svg){html: width=55%}

While the number of subproblems is not limited, usually two or three subproblems suffice. For example, a system of three subproblems could consist of a free flow model at the top coupled to a porous medium in the middle layer and another porous medium with differing equations of state at the bottom. The data exchange would occurr at the coupling boundaries. An example for a system of two subproblems could consist of a free flow layer at the top that is coupled to a pore network model at the bottom. One could also keep the model consistent while varying the dimensionality across the coupling boundary, such as coupling a 3D Darcy model to a 2D Darcy model.

Available implementations include @ref DarcyDarcyCoupling, @ref FreeFlowPoreNetworkCoupling, @ref FreeFlowPorousMediumCoupling and @ref StokesDarcyCoupling.
