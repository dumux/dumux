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
* enalbeThermalDispersion():
    - returns ture, if the model considers thermal dispersion.

# ModelTraits

`ModelTraits` classes are trait classes which encapsulate the main constraints for the model relevant classes. Most important classes affected are @problem,  @assembler and @iofields. The `ModelTraits` define the restrictions to the dervied classes in order to guarantee to some degree that they implemented the required member functions. In general, it is advisable to derive a `ModelTraits` class from a `BasicModelTraits` class.

## Key Members

Note that the following members are inside the namespace `Properties`. The following properties are incorporated in most of the models. To obtain the Property defined in the trait class, you may use `using ExampleProperty = GetPropType<TypeTag, Properties::ExampleProperty>;`. 

* LocalResidual:
    -  local residual is implemented and specifies its type.
* ModelTraits:
    - Specifies the main requirements of the model. If an implementation using this model does not fulfil them, a compiler error should hinder the compilation of the program, to prevent the usage of incorrect models.
    - Note: It is also possible, to use isothermal model traits by extending them with nonisothermal model traits.
* PrimaryVariables:
    - defines the type of `PrimaryVariables` used, most importantly the length of the local solution, which is defined by its submember `PrimaryVariablesVector`. The `PrimaryVariables` class also defines whether primary variable switches are used or not.
* VolumeVariables:
    - trait class specifying requirements for the `VolumeVariables`. Normally, the `VolumeVariables` collect the following properties:
        - `PrimaryVariables` (PV)
        - `FluidSystem` (FSY)
        - `FluidState` (FS)
        - `SolidSystem` (SSY)
        - `SolidState` (SST)
        - `ModelTraits` (MT)
        - `DiscretizationMethod` (DM)
    optionally, there also may be used:
        - `EffectiveDiffusivityModel` (EDM)
        - `MolecularDiffusionType` (DT)
* AdvectionType:
    - Specifies what kind of advection model is implemented, also ensures that an advection model is implemented.
* Formulation:
    - to decide, what kind of formulation is used for the model equations, e.g. pressure saturation.
* FluidState:
    - Ensures that the `FluidState` class is compatible with the `FluidSystem`.
* FluidSystem:
 -  the FluidSystem that is used by the model. The trait class ensures, that the `FluidSystem` can be used for the model.
 * ThermalConductivityModel:
    - for nonisothermal models: Ensures that a model for thermal conductivity is implemented.
