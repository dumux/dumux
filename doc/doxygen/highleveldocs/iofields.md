# IOFields

This class can used to add model specific output data to a custom @ref vtkoutputmodule "vtkWriter". An `IOFields` class can be taylored to a specific `Model` class, see also @ref model.


### Key functionalities

* initOutputModule():
    - Initialize a `vtkOutputModule` by specifying which variables shoud be written out and add the them via the function `addVolumeVariable()` to the `vtkOutputModule`.

### Overview

@mermaid{iofields}