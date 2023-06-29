# Model

`Model` classes are trait classes which are used as basis classes for concrete `Problem` classes. One could say, the `Model` is the *abstraction* of a conceptual model, while the `Problem` is the *concretion* of this conceptual model. The `Model` class defines, what a derived `Problem` class has to incorporate in order to be correct. These constraints should be chosen such, that they ensure derived `Problem` classes do not incorporate contradicting subelements. See also @ref problem.



### Key Members

## Basic Indices

* numEq():
    - the number of equations used in the model.
* numFluidPhases():
    - the number of fluid phases considered in the model.
* enableAdvection():
    define, whether the model considers advection

* enableEnergyBalance():
    define, if the model has an energy balance.

## Specific Indices
* enableMolecularDiffusion():
    for a compositional model: decide if molecular diffusion is modeled.


## Basic Propetries

Note that the following members are inside the namespace `Properties`. The following properties are incorporated in most of the models. To obtain the Property defined in the trait class, you may use `using ExampleProperty = GetPropType<TypeTag, Properties::ExampleProperty>;`. 

* LocalResidual:
    - Ensures that a local residual is implemented and specifies its type.
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
        
## More Specific Properties

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
