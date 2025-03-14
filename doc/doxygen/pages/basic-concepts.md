# Basic concepts

## Properties / Property System
In DuMu<sup>x</sup>, the property system provides a flexible way to configure simulations at compile time. Properties are structs that define types which determine how different parts of the framework work together. The system ensures a consistent class hierarchy by allowing users to inherit from predefined models and customize specific properties as needed.
Users typically collect their property customizations in a `properties.hh` file specific to their simulation setup. For more detailed and technical information about the property system, see @ref Properties.

## Problem
A "problem" in DuMu<sup>x</sup> conceptually characterizes the simulated scenario by specifying initial and boundary conditions, scenario-specific source terms, and model parameters. A problem is a class that defines the problem interface methods the selected model requires. Users implement their own problem class for each simulation scenario. User implementations of the problem typically inherit from one of the basic implementations (@ref Dumux::FVProblem, @ref Dumux::FVProblemWithSpatialParams, @ref Dumux::PorousMediumFlowProblem) and customizing the default implementations by defining specific interface functions like @ref Dumux::FVProblem::source in the derived class. Examples of problem class implementations can be found in the [documented examples](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples), and almost every end-user test application in the `test` folder (typically in a file called `problem.hh`).

### SpatialParams
All problem implementations derived from @ref Dumux::FVProblemWithSpatialParams (for example, @ref Dumux::PorousMediumFlowProblem) contain an instance of a class that specified position-dependent parameters. The type of this class is set by specializing the property `Dumux::Properties::SpatialParams`. Examples of spatial parameter class implementations can also be found in the documented examples and many end-user test applications.

## More basic concepts
This documentation page is still incomplete. Contributing by describing other basic concepts in DuMu<sup>x</sup>.
To go to the source file of this documentation page, follow the edit button below.

