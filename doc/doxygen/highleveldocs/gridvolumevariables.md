@page gridvolumevariables gridVolumeVariables

In general there are two options in Dumux:
1. Caching enabled
2. Caching disabled

If caching is enabled, a gridVolumeVariables object (e.g., curGridVolVars) creates a vector of @ref volumeVariables.
Since the smallest entity of a volume is a sub-control-volume in Dumux, the length of the vector is the number of sub-control-volumes.
One can access entities of that vector with the respective sub-control-volume-index.
The localAssembler gets access to all needed entities of that vector. This includes all sub-control-volumes of that element.
In that case the elementVolumeVariables only forward the entries of the globally stored vector.

In the case of disabled caching, there is no global vector created.

- @subpage elementvolumevariables