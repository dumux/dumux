## LocalResidual {#localresidual}
<!-- @page localResidual LocalResidual -->


The localResidual contains the actual PDE of the problem one wants to solve. In general, one is interested in PDE's of the form

\f{equation}{
 \underbrace{\frac{\partial S(\boldsymbol{u})}{\partial t}}_{1)}
- \underbrace{\nabla \cdot \boldsymbol{F(\boldsymbol{u})} }_{2)}
- \underbrace{q}_{3)}
= 0
\f}
where 1) describes the storage term which handles changes over time, 2) stands for the flux term, which can have advective or diffusive character and 3) is the source/sink term.

The localResidual is defined as LHS-RHS evaluated for one sub-control-volume.

\f{equation}{
r(\boldsymbol{u})
= \underbrace{\frac{\partial S(\boldsymbol{u})}{\partial t}}_{1)}
- \underbrace{\nabla \cdot \boldsymbol{F(\boldsymbol{u})} }_{2)}
- \underbrace{q(\boldsymbol{u})}_{3)}
\f}


### Functionalities

1. evalFluxAndSource()
   1. evalSource()
   2. evalFlux()
2. evalStorage()

### Overview

@mermaid{localresidual}
