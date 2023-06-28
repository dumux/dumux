@page volumevariables VolumeVariables

volumeVariables is a class that gives access and lets you store variables are actually needed for the computation. Here, one could implement a function that calculates the density from primaryVariables.
However, in DuMuX there exists another layer. For instance, access to variables that are connected to your fluid of interest can be accessed and stored via the @subpage fluidstate.
If you keep using the example of variables that are connected to your fluid of interest, equations to calculate the variables are per default defined in the @subpage fluidsystem.
The same logic applies to other scenarios, for instance variables that are connected to the solid (i.e., porosity).
Per default, the VolumeVariables class acts gives you access to the next layer. However, this is not mandatory, one could implement everything here.

### Key functionalites

- completeFluidState()
  - sets the properties of the @ref fluidstate, that are primaryVariables
  - calculates the rest of the variables concerning the fluid using the functions defined in the @ref fluidsystem


### Overview

@mermaid{volumevariables}
