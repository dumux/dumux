@page elementvolumevariables ElementVolumeVariables

In short, elementVolumeVariables are the view of a element towards all @subpage volumevariables located at sub-control-volumes in it's respective stencil.

In the case of enabled caching, the elementVolumeVariables forward the entries belonging to that element from the globally stored vector described in @ref gridvolumevariables

If caching is disabled, the elementVolumeVariables do create a vector of @ref volumeVariables for each element on the fly.

In both cases the elementVolumeVariables will be forwarded to the @ref localassembler. Since the localassembler only needs information of the sub-control-volumes in the respective stencil.

### Key functionalities

- bind()
  - for caching enabled, precomputes dirichlet boundary values in the stencil
  - for caching disabled, computes all @ref volumevariables in the stencil.
