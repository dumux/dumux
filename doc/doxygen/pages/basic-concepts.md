# Basic concepts

## Properties / Property System
In DuMuX, the property system provides a flexible way to configure simulations at compile time. Properties are structs that define types which determine how different parts of the framework work together. The system ensures a consistent class hierarchy by allowing users to inherit from predefined models and customize specific properties as needed.
Users typically collect their property customizations in a `properties.hh` file specific to their simulation setup. For more detailed and technical information about the property system, see @ref Properties.

## SpatialParams
The SpatialParams define parameters which can be dependent on the location in space within the computational domain. For a Darcy-scale porous-medium simulation, these are typically porosity and permeability as well as parameters appearing in the constitutive relations, such as the capillary entry pressure. The specified spatial parameters are usually stored in a `spatialparams.hh` file. More information can be found @ref SpatialParameters.
