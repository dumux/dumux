# IOFields

This class can used to add model specific output data to a custom @ref vtkoutputmodule "vtkWriter". An `IOFields` class can be taylored to a specific `Model` class, see also @ref model. 


### Key functionalities

* initOutputModule():
    - Initialize a `vtkOutputModule` by specifying which variables shoud be written out and add the them via the function `addVolumeVariable()` to the `vtkOutputModule`.
* primaryVariableName(int pvIdx, int state)
    - returns the human readable name for the primary variable with index `pvIdx`. Models that do not use primary variable switching, will have always `state = 0`. In case of a primary variable switch (PVS) model, `state` has to be explicitly given when calling. Note that the state has to be known to correctly interpret the primary variables in case of a PVS model.
### Overview

@mermaid{iofields}
