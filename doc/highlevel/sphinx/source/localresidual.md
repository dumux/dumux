## localresidual

* Is a representation of the discretized differential equations - the differential equations are implented in residual form such as

```{math}
\underbrace{\varPhi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}}_{1}
- \nabla \cdot \left( \varrho_\alpha \frac{k_{r,\alpha}}{\mu_\alpha} \boldsymbol{K}\left( \nabla p_\alpha - \varrho_\alpha \boldsymbol{g} \right) \right)
- q_\alpha
= 0
```

which describes the mass balance of a multi-phase model. Or in general

```{math}
r(\boldsymbol{u})
= \underbrace{\frac{\partial S(\boldsymbol{u})}{\partial t}}_{1)}
- \underbrace{\nabla \cdot \boldsymbol{F(\boldsymbol{u})} }_{2)}
- \underbrace{q}_{3)}
= 0
```

where 1) describes the storage term which handles changes over time, 2) stands for the flux term, which can have advective or diffusive character and 3) is the source/sink term for one specific sub-control volume.

* The individual terms of the local residual are computed for EACH sub-control volume and the values are stored in a vector which has an entry for each sub-control volume.

* Called by: LocalAssembler to compute entries in the jacobian and residual

* Source term is basically a right-hand side, can be used for anything that cannot be described within the storage or flux term


### Functionalities
* evalStorage()
* computeSource()
* computeFlux()

### Connected Classes

```{mermaid}

   flowchart LR
   1[localAssembler] --> 2[localResidual]
   click 1 "./localassembler.html"
```