#!/usr/bin/env python3

#  Compile some C++ code for reference
from dumux.common import *
from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

# debugging/testing code
def PrintProblemTest(problem):
    includes = problem._includes + ["test/python/test_boundaryloop.hh"]
    typeName = "Dumux::Python::PrintProblemTest<{}>".format(problem._typeName)
    moduleName = moduleName = "printbs_" + hashIt(problem._typeName)
    generator = SimpleGenerator("PrintProblemTest", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName)
    return module.PrintProblemTest(problem)

############################################################
# The actual Python test
############################################################
from dune.grid import structuredGrid
from dumux.discretization import GridGeometry
from dumux.common import BoundaryTypes, FVProblem

gridView = structuredGrid([0,0,0],[1,1,1],[3,3,3])

gridGeometry = GridGeometry(gridView, discMethod="box")
gridGeometry.update()

@FVProblem(gridGeometry)
class Problem:
    numEq = 2
    name = "python_problem"

    def boundaryTypes(self, element, scv):
        bTypes = BoundaryTypes(self.numEq)
        bTypes.setDirichlet()
        return bTypes

    def dirichlet(self, element, scv):
        if scv.center()[0] > 0.5:
            return [0.5, 0.5]
        else:
            return [1.0, 0.0]

    def sourceAtPos(self, globalPos):
        return [globalPos[0]]

problem = Problem()
print("Name of the problem: {}".format(problem.name))
print("-- Number of equations: {}".format(problem.numEq))
print()

# test C++/Python compatibility of the hybrid class

# print some info from C++
printer = PrintProblemTest(problem)
printer.print()

# print the same info from Python
numNeumann = 0
numDirichlet = 0
totalSource = 0
for e in gridView.elements:
    fvGeometry = problem.gridGeometry().localView()
    fvGeometry.bind(e)
    for scv in fvGeometry.scvs():
        bTypes = problem.boundaryTypes(element=e, scv=scv)
        if bTypes.isDirichlet():
            numDirichlet += 1
        elif bTypes.isNeumann():
            numNeumann += 1
        totalSource += problem.sourceAtPos(scv.center())[0]*scv.volume()

print("[python] Found {} Neumann faces and {} Dirichlet faces".format(numNeumann, numDirichlet))
print("[python] Total source {:.2f} kg/s".format(totalSource))
