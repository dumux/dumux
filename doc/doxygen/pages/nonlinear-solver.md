# Non-linear solver: Newton
The Newton solver is used to solve a non-linear system of equations.
In DUMuX, users have the flexibility to adjust various Newton solver parameters to optimize convergence and efficiency. These aparmeters can be set in the input file, e.g., `params.input`.

## Parameters

* UseLineSearch
* EnableChop
* EnablePartialReassembly
* EnableAbsoluteResidualCriterion
* EnableShiftCriterion
* EnableResidualCriterion
* SatisfyResidualAndShiftCriterion
* MaxRelativeShift
* MaxAbsoluteResidual
* ResidualReduction
* MinSteps
* MaxSteps
* TargetSteps
* ReassemblyMinThreshold
* ReassemblyMaxThreshold
* ReassemblyShiftWeight
* RetryTimeStepReductionFactor
* MaxTimeStepDivisions
