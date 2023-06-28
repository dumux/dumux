@page localassembler LocalAssembler

The localAssembler uses the corresponding @subpage localresidual to calculate the Jacobian and the residual-vector.
It only works on the element view and hence, only gets acces to the variables in the stencil (e.g., @ref elementvolumevariables).
If caching is disabled, the variables in the stencil first need to be build. Since the @ref localresidual works on a local view, the LocalAssembler also needs to map the @ref localresidual from the local scv indices to the global scv indices.

### Key functionalities

- evalFluxAndSourceResidual()
  - calls evalFluxAndSource() of the @ref localresidual.
- evalStorageResidual()
  - calls evalStorage() of the @ref localresidual.

### Overview

@mermaid{localassembler}
