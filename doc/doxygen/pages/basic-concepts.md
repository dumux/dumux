# Basic concepts

## Properties / Property System
In DuMuX, the property system provides a flexible way to configure simulations at compile time. Properties are structs that define types which determine how different parts of the framework work together. The system ensures a consistent class hierarchy by allowing users to inherit from predefined models and customize specific properties as needed.
Users typically collect their property customizations in a `properties.hh` file specific to their simulation setup. For more detailed and technical information about the property system, see @ref Properties.

## Problem
A Problem in DuMux represents the conceptual framework where the scenario being simulated is characterized through the specification of initial and boundary conditions, as well as source terms. It is recommended to save all specifications of the initial and boundary conditions in one file e.g. `problem.hh`.
