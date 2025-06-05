@addtogroup MultiDomain

The MultiDomain framework is an abstract concept that allows coupling different pyhsical models and/or computational domains to each other. Here, we will only briefly present the concept while more details can be found in \cite Koch2021.

![Interaction regions between different subproblems.](multidomainSubproblemsInteractions.svg){html: width=40%}

In DuMu<sup>x</sup>, each subproblem embodies the conceptual representation of a physical model applied to some computational domain. In a non-coupled context, simulations usually consist of exactly one problem entity. Here, for illustration purposes, we will assume to have three distinct subproblems SP1, SP2 and SP3. The MultiDomain framework allows to have several subproblems which can also interact with each other. Parts of the computational domain might be isolated to a single subproblem while other regions might be shared by multiple subproblems. The orange, green and purple areas illustrate regions with coupling.

![Global matrix structure of coupled subproblems.](multidomainAssembler.svg){html: width=55%}

The subproblems on their own are conceptually isolated, meaning they usually do not know of the other subproblems. However, for the multi domain approach, DuMu<sup>x</sup> provides a coupling manager class. This class manages the data exchange between subproblems and ensures that coupling terms are communicated across subproblems. Closely tied to the coupling manager, there is also an optional coupling mapper class which may take the responsibility of specifying the so-called coupling stencils to alleviate the struture of the coupling manager. A coupling stencil defines all the grid elements that participate in the computation of the residual for one specific element. During each Newton iteration, the coupling manager communicates with the multidomain assembler to construct the global matrix. For each subproblem, the assembly results in a block matrix on the diagonal of the global matrix which represents the isolated subproblems. The submatrices on the off-diagonal entries represent the interactive regions in which two subproblems are coupled to each other. In case only two subproblems were available, the global matrix would only have two block matrices on the diagonal and two off-diagonal submatrices in total.

![Examples for possible coupled systems.](coupling_examples.svg){html: width=100%}

The MultiDomain abstraction enables developers to couple arbitrary models and computational domains. Currently, available implementations in DuMu<sup>x</sup> encompass the @ref BoundaryCoupling, @ref DualNetworkCoupling, @ref EmbeddedCoupling and finally the @ref FacetCoupling.
