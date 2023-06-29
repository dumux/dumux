# Assembler


The `Assembler` is responsible for calculating the global residual vector and the global Jacobian matrix. It relies on a discretization-specific local assembler engine. The Jacobian matrix is the partial derivative of the residual with respect to each entry of the solution vector. Dumux uses an element-wise assembly algorithm. For each element, a local assembler is instantiated. Typically, the `assembleJacobianAndResidual` method is called which assembles both residual and Jacobian in one go.

### Key functionalites

- assembleJacobianAndResidual()
   - Calls the assembleJacobianAndResidual() function of the `LocalAssembler` for each element.

## LocalAssembler

The local `Assembler` employs the relevant `LocalResidual` to compute the Jacobian and the residual vector. It is limited to the element view and therefore only has access to the variables in the respective stencil, such as element volume variables. If caching is disabled, the `Assembler` must first build the stencil variables. As the `LocalResidual` operates on a local view, the `LocalAssembler` must also map the `LocalResidual` from local `Sub-control-volume`-indices to global `Sub-control-volume`-indices.

### Key functionalities

- evalFluxAndSourceResidual()
  - calls evalFluxAndSource() of the `LocalResidual`.
- evalStorageResidual()
  - calls evalStorage() of the `LocalResidual`.

## LocalResidual

The `LocalResidual` contains the actual PDE of the problem one wants to solve. In general, one is interested in PDE's of the form

```math
 \underbrace{\frac{\partial S(\boldsymbol{u})}{\partial t}}_{1)}
- \underbrace{\nabla \cdot \boldsymbol{F(\boldsymbol{u})} }_{2)}
- \underbrace{q}_{3)}
= 0
```
where 1) describes the storage term which handles changes over time, 2) stands for the flux term, which can have advective or diffusive parts and 3) is the source/sink term.

The `LocalResidual` is defined as LHS-RHS evaluated for one sub-control-volume.

```math
r(\boldsymbol{u})
= \underbrace{\frac{\partial S(\boldsymbol{u})}{\partial t}}_{1)}
- \underbrace{\nabla \cdot \boldsymbol{F(\boldsymbol{u})} }_{2)}
- \underbrace{q(\boldsymbol{u})}_{3)}
```

### Key Functionalities

- evalFluxAndSource()
   - evalSource()
   - evalFlux()
- evalStorage()