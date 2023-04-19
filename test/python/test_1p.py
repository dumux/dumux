#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


import sys

from dune.grid import structuredGrid
from dune.istl import blockVector, CGSolver, SeqSSOR

from dumux.common import initParameters, printParameters, getParam
from dumux.common import BoundaryTypes, Model, Property
from dumux.discretization import GridGeometry, GridVariables
from dumux.assembly import FVAssembler
from dumux.porousmediumflow import (
    PorousMediumFlowProblem,
    PorousMediumFlowVelocityOutput,
    FVSpatialParamsOneP,
)
from dumux.material import FluidSystem, Component
from dumux.io import VtkOutputModule

# Initialize parameters
initParameters(
    argv=sys.argv,
    params={
        "Problem.EnableGravity": True,
        "Vtk.AddVelocity": False,
        "Assembly.NumericDifference.PriVarMagnitude": 1e5,
    },
)

discMethod = getParam("DiscMethod", default="box")
if discMethod not in ["box", "cctpfa"]:
    raise NotImplementedError(discMethod + " not supported yet. Use cctpfa or box")

diffMethod = getParam("DiffMethod", default="numeric")
if diffMethod not in ["analytic", "numeric"]:
    raise NotImplementedError(diffMethod + " must be analytic or numeric")


# Set up the grid and the grid geometry
gridView = structuredGrid([0, 0], [1, 1], [10, 10])
gridGeometry = GridGeometry(gridView=gridView, discMethod=discMethod)

# Set up the model
model = Model(inheritsFrom=["OneP"], gridGeometry=gridGeometry)

# Tell model to use a particular local residual type
model["LocalResidual"] = Property.fromCppType(
    "OnePIncompressibleLocalResidual<TypeTag>",
    cppIncludes=["dumux/porousmediumflow/1p/incompressiblelocalresidual.hh"],
)

# Setup the fluid system
h20 = Component("SimpleH2O")
onePLiquid = FluidSystem("OnePLiquid", component=h20, scalar=model["Scalar"])
model["FluidSystem"] = Property.fromInstance(onePLiquid)


# Define the spatial parameters
@FVSpatialParamsOneP(gridGeometry=gridGeometry)
class SpatialParams:
    dimWorld = gridGeometry.gridView.dimWorld
    lensLowerLeft = [0.2, 0.2]
    lensUpperRight = [0.8, 0.8]

    def isLens(self, globalPos):
        eps = 1.5e-7
        for i in range(self.dimWorld):
            if (globalPos[i] < self.lensLowerLeft[i] + eps) or (
                globalPos[i] > self.lensUpperRight[i] - eps
            ):
                return False
        return True

    def permeability(self, element, scv, elemSol):
        globalPos = scv.dofPosition
        # permeability can be either given
        # as scalar or tensor
        if self.isLens(globalPos):
            return [[1e-12, 0], [0, 1e-12]]

        return 1e-10

    def porosityAtPos(self, globalPos):
        return 0.4


spatialParams = SpatialParams()
model["SpatialParams"] = Property.fromInstance(spatialParams)

# Define the problem
@PorousMediumFlowProblem(gridGeometry=gridGeometry, spatialParams=spatialParams)
class Problem:
    numEq = 1

    def boundaryTypesAtPos(self, globalPos):
        bTypes = BoundaryTypes(self.numEq)
        eps = 1e-6

        if globalPos[1] < eps or globalPos[1] > gridGeometry.bBoxMax[1] - eps:
            bTypes.setDirichlet()
        else:
            bTypes.setNeumann()

        return bTypes

    def dirichletAtPos(self, globalPos):
        dp_dy_ = -1.0e5
        return 1.0e5 + dp_dy_ * (globalPos[1] - gridGeometry.bBoxMax[1])

    def sourceAtPos(self, globalPos):
        return 0.0

    def name(self):
        return "python_problem"

    def paramGroup(self):
        return ""

    def neumann(self, element, fvGeometry, scvf):
        return 0


# and set it as a model property
problem = Problem()
model["Problem"] = Property.fromInstance(problem)

# Initialize the GridVariables and the Assembler
sol = blockVector(gridGeometry.numDofs)
gridVars = GridVariables(problem=problem, model=model)
gridVars.init(sol)
assembler = FVAssembler(problem=problem, gridVariables=gridVars, model=model, diffMethod=diffMethod)
print("num dofs: ", assembler.numDofs)

# Assemble the Jacobian and the residual
assembler.assembleJacobianAndResidual(sol)
print("Assembly done")
res = assembler.residual
jac = assembler.jacobian

# Solve the linear system
S = CGSolver(jac.asLinearOperator(), SeqSSOR(jac), 1e-10)
res *= -1
_, _, converged, _, _ = S(sol, res)

if not converged:
    raise Exception("CGSolver has not converged")

# Write to vtk
testName = "test_1p_" + discMethod + "_" + diffMethod
output = VtkOutputModule(gridVariables=gridVars, solutionVector=sol, name=testName)
velocityoutput = PorousMediumFlowVelocityOutput(gridVariables=gridVars)
output.addVelocityOutput(velocityoutput)
output.addVolumeVariable(lambda vv: vv.pressure(), "p")
output.write(1.0)

# Print used and unused parameters
printParameters()
