@page elementvolumevariables elementVolumeVariables

In the case of enabled caching, the elementVolumeVariables forward the entries belonging to that element from the globally stored vector described in @ref gridvolumevariables

If caching is disabled, the elementVolumeVariables do create a vector of volumeVariables for each element.

In both cases the elementVolumeVariables will be forwarded to the @ref localresidual.

- @subpage volumevariables