from dune.grid import structuredGrid
from dumux.discretization import GridGeometry

dim = 2
cells = 25

gridView = structuredGrid([0]*dim, [1]*dim, [cells]*dim)

gridGeometry = GridGeometry(gridView, discMethod="cctpfa")
gridGeometry.update()

from dumux.common import FVProblem
@FVProblem(gridGeometry)
class Problem:
    numEq = 1
    name = "finitevolume"

    def boundaryTypes(self, element, scv):
        bTypes = BoundaryTypes(self.numEq)
        bTypes.setDirichlet()
        return bTypes

    def dirichlet(self, element, scvf):
        return 0.0

    def intial(self, entity):
        if entity.geometry.center.two_norm < 0.25:
            return 1.0
        else:
            return 0.0
problem = Problem()

from dune.common import FieldVector
v = FieldVector([1, 1])
def f(u):
    return v * u

def g(uIn, uOut, n):
    if (v * n) > 0:
        fu = f(uIn)
    else:
        fu = f(uOut)
    return (fu * n) * scvf.area()

import numpy as np
uh = np.zeros(gridGeometry.numDofs())

from dune.grid import gridFunction
@gridFunction(gridView)
def solution(element, x):
    fvGeometry = gridGeometry.localView()
    fvGeometry.bind(element)
    for scv in fvGeometry.scvs():
        return uh[scv.dofIndex()]

def plot():
  solution.plot(clim=[0,1+1e-6], gridLines=None)

# Initialize uh
for e in gridView.elements:
    fvGeometry = gridGeometry.localView()
    fvGeometry.bind(e)
    for scv in fvGeometry.scvs():
        uh[scv.dofIndex()] = problem.initial(e)
plot()

from dumux.common import TimeLoop
dt = 0.5 / (cells * v.two_norm) # CFL condition
timeloop = TimeLoop(startTime=0.0, dt=dt, endTime=1.0, verbose=True)
timeloop.setPeriodicCheckPoint(0.25)

timeloop.start()
while not timeloop.finished():
    upd = np.zeros(gridGeometry.numDofs())

    for e in gridView.elements:
        fvGeometry = gridGeometry.localView()
        fvGeometry.bind(e)

        for scvf in fvGeometry.scvfs():
            i = scvf.insideScvIdx()
            uIn = uh[i]

            if scvf.boundary:
                uOut = problem.dirichlet(e, scvf)[0]
            else:
                j = scvf.outsideScvIdx()
                uOut = uh[j]

            flux = g(uIn, uOut, scvf.unitOuterNormal())
            upd[i] -= flux

            if not scvf.boundary():
                j = scvf.outsideScvIdx()
                upd[j] += flux

        for scv in fvGeometry.scvs():
            upd[scv.dofIndex()] /= scv.volume()

    uh += dt * upd
    timeloop.advanceTimeStep()
    timeloop.reportTimeStep()

    if timeloop.isCheckPoint():
      plot()

timeloop.finalize()
