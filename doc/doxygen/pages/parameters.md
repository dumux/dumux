# Runtime Parameters

Simulation parameters can be parsed to the program via a parameter file or via the command line.

After having run the example application from the getting started guide you will get the following output at the end of the simulation run. If you did not get the output, call Dumux::Parameters::print() in your main file.

```bash
# Runtime-specified parameters used:
[ Grid ]
Cells = "48 32"
UpperRight = "6 4"
[ Newton ]
EnablePartialReassembly = "true"
[ Problem ]
EnableGravity = "true"
Name = "2p"
[ SpatialParams ]
LensLowerLeft = "1.0 2.0"
LensUpperRight = "4.0 3.0"
[ TimeLoop ]
DtInitial = "250"
TEnd = "3000"

# Global default parameters used:
[ Assembly ]
NumericDifferenceMethod = "1"
[ Flux ]
UpwindWeight = "1.0"
[ LinearSolver ]
MaxIterations = "250"
ResidualReduction = "1e-13"
Verbosity = "0"
[ LinarSolver.Preconditioner ]
Iterations = "1"
Relaxation = "1.0"
[ Newton ]
EnableAbsoluteResidualCriterion = "false"
EnableChop = "false"
EnableResidualCriterion = "false"
EnableShiftCriterion = "true"
MaxAbsoluteResidual = "1e-5"
MaxRelativeShift = "1e-8"
MaxSteps = "18"
MinSteps = "2"
ResidualReduction = "1e-5"
SatisfyResidualAndShiftCriterion = "false"
TargetSteps = "10"
UseLineSearch = "false"
[ TimeLoop ]
MaxTimeStepSize = "1e300"
[ Vtk ]
AddProcessRank = "true"
AddVelocity = "false"

# Unused parameters:
Grid.LowerLeft = "0 0"
```

A number of things can be learned:

* run-time parameters can be changed without re-compiling
* default parameters are set by default
* unused parameters are not used by the simulation (maybe typo or wrong group in input file)

## Parameter Values
To get the value of an input parameter please use:

```c++
static const TYPE paramname = getParam<TYPE>("GROUPNAME.PARAMNAME");
```

If you also want to set a default value for a parameter, just add it like this:

```c++
static const TYPE paramname = getParam<TYPE>("GROUPNAME.PARAMNAME", default);
```

As this function call is relatively expensive, the respective variables should always be `static` (e.g., if used in a loop). When dealing with multiple group names, e.g., in the context of coupled models, the following methods might be more convenient:

```c++
auto modelParamGroup0 = "Model0";
static const TYPE paramname0 = getParamFromGroup<TYPE>(modelParamGroup0, "GROUPNAME.PARAMNAME");
auto modelParamGroup1 = "Model1";
static const TYPE paramname1 = getParamFromGroup<TYPE>(modelParamGroup1, "GROUPNAME.PARAMNAME");
```

The Dumux::FVProblem class provides a convenience function Dumux::FVProblem::paramGroup.

The parameters can then be specified in the input file:

```ini
[ Model0.Grid ]
File = file0.dgf
[ Model1.Grid ]
File = file1.dgf
```
