# VtkOutputModule

Strictly, a VtkOutputModule is not necessary to run a simulation in DuMu<sup>x</sup>. However, is used for writing DuMu<sup>x</sup> simulation data to VTK format, which can then be plotted via a software e.g. [ParaView](https://www.paraview.org/). All VtkOutputModules inherit from the [base VtkOutputModule](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/io/vtkoutputmodule.hh). When adding custom data fields for plotting, it is important to keep in mind that these data fields cannot be temporary containers, i.e. the data that is passed to the writer needs to outlive the instance of the writer. One can either choose to instantiate a `vtkWriter` from the class `VtkOutputModuleBase` or from a predefined class `VtkOutputModule`. Using the base module also requires adding model specific output fields via the class @subpage iofields "IOFields", which adds custom data fields via the function `addVolumeVariable()`.


### Key functionalities

* addField(vector, name, ...):
    - Add a scalar or vector valued vtk field `vector` to the tracked output quantities of the `vtkOutputModule`.
* write(time):
    - Write the data for this timestep to file in four steps: (1) assemble all registered variable fields, (2) register them with the vtk writer, (3) the writer writes the output into files, (4) clear the writer for the next time step
* addVelocityOutput():
    - Add a velocity output policy. For example this toggles whether the darcy velocity for porous medium models should be included in the output or not.
* addVolumeVariable(variable):
    - Similar to `addField()`, but instead of adding a vtk field, this function adds a scalar or vector valued object of @ref volumevariables.

### Overview

@mermaid{vtkoutputmodule}
