#!/usr/bin/env python3

import dune
from dune.grid import structuredGrid, gridFunction, OutputType, P1VTKFunction
from dune.istl import blockVector, BlockVector, CGSolver, SeqJacobi

from dumux.discretization import GridGeometry, GridVariables
from dumux.assembly import FVAssembler
from dumux.common.properties import TypeTag
from dumux.common import BoundaryTypes, FVProblem, Parameters
from dumux.material.fluidsystems import FluidSystem
from dumux.material.components import Component, listComponents
from dumux.material.spatialparams import SpatialParams

# Initialize the paramaters
parameters = Parameters()
parameters.init('params.input')

# Choose the discretization method
discMethod = 'box'
# discMethod = "cctpfa"

# Set up the grid and the grid geometry
gridView = structuredGrid([0,0],[1,1],[15,15])
gridGeometry = GridGeometry(gridView, discMethod=discMethod)
gridGeometry.update()

# this TypeTag is for testing purposes only
testTypeTag = TypeTag("Test")
testTypeTag["UseMoles"] = False

# our model TypeTag
myModel = TypeTag('MyModel', inheritsFrom=[testTypeTag,'OneP', ('BoxModel' if discMethod == "box" else 'CCTpfaModel')])

# set some other TypeTag
class MyScalar:
        def __init__(self):
            self._typeName = 'double'

myModel['Scalar'] = MyScalar()

class MySpecialProperty:
    def __init__(self):
        self._typeName = 'Dumux::CubicSpline<Scalar>'
        self._includes = ['<dumux/common/cubicspline.hh>']
        self._requiredPropertyTypes = ['Scalar']

class MyGrid:
    def __init__(self):
        self._typeName = 'typename ' + gridView._typeName + '::GridView::Grid'
        self._includes = ['dune/grid/yaspgrid.hh']

myModel['Grid'] = MyGrid()

class MyLocalResidual:
    def __init__(self):
        self._typeName = 'OnePIncompressibleLocalResidual<TypeTag>'
        self._includes = ['dumux/porousmediumflow/1p/incompressiblelocalresidual.hh']

myModel['LocalResidual'] = MyLocalResidual()

spatialParams = SpatialParams(gridGeometry, MyScalar())
myModel['SpatialParams'] = spatialParams

h20 = Component("SimpleH2O")
onePLiquid = FluidSystem(MyScalar(), h20)
myModel['FluidSystem'] = onePLiquid

# define the Problem
@FVProblem(gridGeometry, spatialParams)
class Problem:
    numEq = 1
    name = "python_problem"

    def boundaryTypes(self, element, scv):
        bTypes = BoundaryTypes(self.numEq)
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
        if discMethod == 'box' and scv.dofIndex() == 134:
            return 1e-1

        if discMethod == 'cctpfa' and scv.dofIndex() == 50:
            return 1e-1

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
myModel['Problem'] = problem

# print the properties
print(myModel.getProperties())

# initialize the GridVariables and the Assembler
gridVars = GridVariables(problem, myModel)
diffMethod = 'analytic'
assembler = FVAssembler(problem, gridVars, myModel, diffMethod)
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

# Write to vtk
@gridFunction(gridView)
def elementGridFunction(element, x):
    elementIdx = gridView.indexSet.index(element)
    return sol[elementIdx]

@gridFunction(gridView)
def vertexGridFunction(element, x):
    if element.type == dune.geometry.cube(2):
        bary = (1-x[0])*(1-x[1]), x[0]*(1-x[1]), (1-x[0])*x[1], x[0]*x[1]
    elif element.type == dune.geometry.simplex(2):
        bary = 1-x[0]-x[1], x[0], x[1]
    idx = gridView.indexSet.subIndices(element, 2)
    return sum(b * sol[i] for b, i in zip(bary, idx))

if discMethod == 'box':
    gridView.writeVTK(problem.name + '_box',
                      pointdata={"solution": vertexGridFunction},
                      outputType=OutputType.ascii)
else:
    gridView.writeVTK(problem.name + '_cctpfa',
                      celldata={"solution": elementGridFunction},
                      outputType=OutputType.ascii)

listComponents()
