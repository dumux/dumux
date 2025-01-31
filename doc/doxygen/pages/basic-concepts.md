# Basic concepts

## Properties/ Property System

In DuMuX, the property system provides a flexible way to configure simulations at compile time.
Properties are structs that define types that determine how different parts of the framework work together. Users can customize their simulations by inheriting from predefined model and modifying specific properties. By inheriting from a predefined model, a consistent tree a classtypes is ensured.
A customization of the predefined properties typically collected in a properties.hh file, specific to a simulation set. 
For more detailed and technical information about the property system @ref Properties .