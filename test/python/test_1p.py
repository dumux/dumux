#!/usr/bin/env python3

import sys

from dune.grid import structuredGrid
from dune.istl import blockVector, CGSolver, SeqSSOR

from dumux.discretization import GridGeometry, GridVariables
from dumux.assembly import FVAssembler
from dumux.common import BoundaryTypes, Parameters, Model, Property
from dumux.porousmediumflow import PorousMediumFlowProblem, PorousMediumFlowVelocityOutput
from dumux.material import FluidSystem, Component, OnePSpatialParams
from dumux.io import VtkOutputModule

try:
    discMethod = sys.argv[1]
    diffMethod = sys.argv[2]
except IndexError:
    print("No discretization method and differentiation method given. Defaulting to box and numeric diff.")
    discMethod = "box"
    diffMethod = "numeric"

if discMethod not in ["box", "cctpfa"]:
    raise NotImplementedError(discMethod + " not supported yet. Use cctpfa or box")

if diffMethod not in ["analytic", "numeric"]:
    raise NotImplementedError(diffMethod + " must be analytic or numeric")

# Initialize the paramaters
parameters = Parameters({
    "Problem.EnableGravity": True,
    "SpatialParams.Porosity": 0.3,
    "SpatialParams.Permeability": 1e-8,
    "Vtk.AddVelocity": False,
    "Assembly.NumericDifference.PriVarMagnitude": 1e5,
})

# Set up the grid and the grid geometry
gridView = structuredGrid([0,0], [1,1], [10,10])
gridGeometry = GridGeometry(gridView=gridView, discMethod=discMethod)

# Set up the model
model = Model(inheritsFrom=['OneP'], gridGeometry=gridGeometry)

# Tell Dumux to use a particular local residual type
model['LocalResidual'] = Property(
    type='OnePIncompressibleLocalResidual<TypeTag>',
    includes=['dumux/porousmediumflow/1p/incompressiblelocalresidual.hh']
)

@OnePSpatialParams(gridGeometry=gridGeometry)
class SpatialParams:
    dimWorld = gridGeometry.gridView.dimWorld
    lensLowerLeft = [0.2, 0.2]
    lensUpperRight = [0.8, 0.8]

    def isLens(self, globalPos):
        eps = 1.5e-7
        for i in range(self.dimWorld):
            if (globalPos[i] < self.lensLowerLeft[i] + eps) or (globalPos[i] > self.lensUpperRight[i] - eps):
                return False
        return True

    def permeability(self, element, scv, elemSol):
        globalPos = scv.dofPosition
        # permeability can be either given
        # as scalar or tensorial value
        if self.isLens(globalPos):
            return [[1e-12, 0],
                    [0, 1e-12]]
        else:
            return 1e-10

    def porosityAtPos(self, globalPos):
        return 0.4

spatialParams = SpatialParams()
model['SpatialParams'] = Property(object=spatialParams)

h20 = Component(name="SimpleH2O")
onePLiquid = FluidSystem(type="OnePLiquid", component=h20, scalar=model['Scalar'])
model['FluidSystem'] = Property(object=onePLiquid)

# define the Problem
@PorousMediumFlowProblem(gridGeometry, spatialParams)
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
        dp_dy_ = -1.0e+5
        return  1.0e+5 + dp_dy_*(globalPos[1] - gridGeometry.bBoxMax[1])

    def sourceAtPos(self, globalPos):
        return 0.0

    def extrusionFactor(self, element, scv):
        return 1.0

    def temperatureAtPos(self, globalPos):
        return 283.15

    def name(self):
        return "python_problem"

    def paramGroup(self):
        return ""

    def neumann(self, element, fvGeometry, scvf):
        return 0

problem = Problem()
model['Problem'] = Property(object=problem)

# initialize the GridVariables and the Assembler
gridVars = GridVariables(problem=problem, model=model)
assembler = FVAssembler(problem=problem, gridVariables=gridVars, model=model, diffMethod=diffMethod)
sol = blockVector(assembler.numDofs)
gridVars.init(sol)
assembler.updateGridVariables(sol)
print("numdofs", assembler.numDofs)

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
output.addVolumeVariable(lambda vv : vv.pressure(), "p")
output.write(1.0)
