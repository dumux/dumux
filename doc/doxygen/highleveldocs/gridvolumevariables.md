@page gridvolumevariables GridVolumeVariables

As the name suggest, the gridVolumeVariables refer to variables that are located inside volumes on your grid.

In general, there are two options in Dumux:
1. Caching enabled
2. Caching disabled

If caching is enabled, a gridVolumeVariables instance (e.g., curGridVolVars) has a vector of @ref volumevariables.
Since the smallest entity of a volume is a sub-control-volume in Dumux, the length of the vector is the number of sub-control-volumes.
One can access entities of that vector with the respective sub-control-volume-index.

If caching is disabled, one can acces the @ref volumevariables only via the @subpage elementvolumevariables.

Exemplary for the two implementations (e.g., caching enabled and disabled) in DuMuX, the implementation of @ref gridvolumevariables is sketched:
If caching is enabled, this is exactly how gridVolumeVariables are implemented. If caching is disabled, the global vector on the left does not exist.


![](gridVariables.svg)

### Key functionalities

If caching is enabled:
<ul>
  <li>update()</li>
  <ul style="list-style:none">
    <li>for every scv you will update the volumeVariables</li>
  </ul>
  <li>volVars()</li>
  <ul style="list-style:none">
    <li>gives you acces to volumeVariables stored at specific entry of the global gridVolumeVariables vector.
    </li>
  </ul>
</ul>

### Overview

@mermaid{gridvolumevariables}
