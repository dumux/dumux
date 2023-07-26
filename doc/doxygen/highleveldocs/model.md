# Model

# BaseModelTraits

The `BaseModelTraits` define the most basic requirements a model has to fulfil, normally this would be the number of model equations, the number of fluid phases, the number of chemical components etc.

### `BasicModelTraits` Key Members
* numEq():
    - returns the number of equations used in the model.
* numFluidPhases():
    - returns the number of fluid phases used in the model.
* numFluidComponents():
    - returns the number of fluid components used in the model.
* numSolidComponents():
    - returns the number of solid components used in the model.
* enableAdvection():
    - returns true if the model considers advection.
* enableMolecularDiffusion():
    - returns true if the model considers molecular diffusion.
* enableEnergyBalance():
    - returns true, if the model incorporates an energy balance. 
* enableThermalDispersion():
    - returns ture, if the model considers thermal dispersion.

# ModelTraits

`ModelTraits` classes are trait classes which encapsulate the main constraints for the model relevant classes. Most important classes affected are @problem,  @assembler and @iofields. The `ModelTraits` define the restrictions to the dervied classes in order to guarantee to some degree that they implemented the required member functions. In general, it is advisable to derive a `ModelTraits` class from a `BasicModelTraits` class.

