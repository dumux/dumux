#!/usr/bin/env python3

import dune
from dune.grid import structuredGrid, gridFunction, OutputType, P1VTKFunction
from dune.istl import blockVector, BlockVector, CGSolver, SeqJacobi

from dumux.discretization import GridGeometry, GridVariables
from dumux.assembly import FVAssembler
from dumux.common.properties import TypeTag, Property
from dumux.common import BoundaryTypes, FVProblem, Parameters
from dumux.material.fluidsystems import FluidSystem
from dumux.material.components import Component, listComponents
from dumux.material.spatialparams import SpatialParams
from dumux.io import VtkOutputModule

# Initialize the paramaters
parameters = Parameters()

paramsDict = {"Problem.EnableGravity":"true"}
# parameters.init('params.input', paramsDict)
parameters.init('params.input', {"Problem.EnableGravity":"true"})

# Choose the discretization method
discMethod = 'box'
# discMethod = "cctpfa"

# Set up the grid and the grid geometry
gridView = structuredGrid([0,0],[1,1],[15,15])
gridGeometry = GridGeometry(gridView, discMethod=discMethod)
gridGeometry.update()

# this TypeTag is for testing purposes only
testTypeTag = TypeTag("Test")
testTypeTag["UseMoles"] = Property(value=False)

# TODO automatically forward declare unknown properties
# testTypeTag["MyScalarValueProp"] = Property(value=123.0)
# testTypeTag["MyIntValueProp"] = Property(value=123)

# our model TypeTag
myModel = TypeTag('MyModel', inheritsFrom=[testTypeTag,'OneP', ('BoxModel' if discMethod == "box" else 'CCTpfaModel')])

myModel['Scalar'] = Property(type='double')
myModel['Grid'] = Property(type='typename ' + gridView._typeName + '::GridView::Grid', includes=['dune/grid/yaspgrid.hh'])
myModel['LocalResidual'] = Property(type='OnePIncompressibleLocalResidual<TypeTag>', includes=['dumux/porousmediumflow/1p/incompressiblelocalresidual.hh'])

spatialParams = SpatialParams(gridGeometry, myModel['Scalar'])
myModel['SpatialParams'] = Property(object=spatialParams)

h20 = Component("SimpleH2O")
onePLiquid = FluidSystem("OnePLiquid", h20, myModel['Scalar'])
myModel['FluidSystem'] = Property(object=onePLiquid)

# define the Problem
@FVProblem(gridGeometry, spatialParams)
class Problem:
    numEq = 1
    name = "python_problem"

    def boundaryTypes(self, element, scv):
        bTypes = BoundaryTypes(self.numEq)
        bTypes.setNeumann()
        if scv.dofPosition()[1] > 1 -1e-8:
            bTypes.setDirichlet()
        return bTypes

    def dirichlet(self, element, entity):
        return 0.0

    def sourceAtPos(self, globalPos):
        return 0.0

    def extrusionFactor(self, element, scv):
        return 1.0

    def temperatureAtPos(self, globalPos):
        return 300.0

    def source(self, element, fvGeometry, scv):
        # if discMethod == 'box' and scv.dofIndex() == 134:
        #     return 1e-1

        # if discMethod == 'cctpfa' and scv.dofIndex() == 50:
        #     return 1e-1

        return 0

    def scvPointSources(self, element, fvGeometry, scv):
        return 0.0

    def paramGroup(self):
        return ""

    def neumann(self, element, fvGeometry, scvf):
        return 0

    def addSourceDerivatives(self, block, element, fvGeometry, scv):
        return 0

problem = Problem()
print("Name of the problem: {}".format(problem.name))
print("Includes of the problem: {}".format(problem._includes))
myModel['Problem'] = Property(object=problem)

# print the properties
print(myModel.getProperties())

# initialize the GridVariables and the Assembler
gridVars = GridVariables(problem, myModel)
assembler = FVAssembler(problem, gridVars, myModel, diffMethod='analytic')
sol = blockVector(assembler.numDofs())
gridVars.init(sol)
assembler.updateGridVariables(sol)
print("numdofs", assembler.numDofs())

# Assemble the Jacobian and the residual
assembler.assembleJacobianAndResidual(sol)
res = assembler.residual()
jac = assembler.jacobian()

# Solve the linear system
S = CGSolver(jac.asLinearOperator(), SeqJacobi(jac), 1e-10)
res *= -1
_, _, converged, _, _ = S(sol, res)

if not converged:
    raise Exception("CGSolver has not converged")

listComponents()

def pressure(volVars):
    return volVars.pressure()

# Write to vtk
output = VtkOutputModule(myModel, gridVars, sol, "test")
output.addField(sol, "x")
output.addVolumeVariable(pressure, "p")
output.addVolumeVariable(lambda vv : vv.density(), "rho")
output.addVolumeVariable(lambda vv : vv.temperature(), "T")
output.addVolumeVariable(lambda vv : vv.saturation(), "S")
output.write(1.0)
