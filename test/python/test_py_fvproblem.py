#!/usr/bin/env python3

import time
from dune.grid import structuredGrid
from dumux.discretization import GridGeometry
from dumux.common import BoundaryTypes, FVProblem, PrintBoundaryStuff

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

problem = Problem()
print("Name of the problem: {}".format(problem.name))
print("-- Number of equations: {}".format(problem.numEq))
print()

# test C++/Python compatibility of the hybrid class

# print some info from C++
printer = PrintBoundaryStuff(problem)
printer.print()

# print the same info from Python
start = time.time()
for e in gridView.elements:
    fvGeometry = problem.gridGeometry().localView()
    fvGeometry.bind(e)
    for scv in fvGeometry.scvs():
        bTypes = problem.boundaryTypes(element=e, scv=scv)
        print("-- scv at {}: isNeumann: {}, isDirichlet: {}".format(scv.dofPosition(), bTypes.isNeumann(), bTypes.isDirichlet()))
        if bTypes.isDirichlet():
            print("  -- Dirichlet values: {}".format(problem.dirichlet(element=e, scv=scv)))
