#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


from dune.common import FieldVector
from dune.grid import structuredGrid, gridFunction, OutputType
from dumux.common import FVProblem, BoundaryTypes, TimeLoop
from dumux.discretization import GridGeometry
import numpy as np

plotting = True
try:
    import matplotlib.pyplot as plt

    plt.ion()
except ImportError:
    print("Warning: Plots are not generated as matplotlib could not be found.")
    plotting = False

########################
# Create grid geometry #
########################

dimension = 2
cells = 20

gridView = structuredGrid([0] * dimension, [1] * dimension, [cells] * dimension)

gridGeometry = GridGeometry(gridView, discMethod="cctpfa")

elementMapper = gridView.indexSet


################################################
# Define problem (initial/boundary conditions) #
################################################


@FVProblem(gridGeometry=gridGeometry)
class Problem:
    numEq = 1

    def name(self):
        return "finitevolume"

    def boundaryTypes(self, element, scv):
        bTypes = BoundaryTypes(self.numEq)
        bTypes.setDirichlet()
        return bTypes

    def dirichlet(self, element, scvf):
        if scvf.ipGlobal.two_norm < 0.25:
            return 1.0
        else:
            return 0.0

    def initial(self, entity):
        return 0.0


problem = Problem()

######################
# Transport equation #
######################

velocity = FieldVector([1] * dimension)
upwindWeight = 1.0


def advectiveFlux(insideConcentration, outsideConcentration, normal):
    normalVelocity = velocity * normal
    upwindConcentration = insideConcentration
    downwindConcentration = outsideConcentration
    if normalVelocity < 0.0:
        upwindConcentration, downwindConcentration = downwindConcentration, upwindConcentration
    return normalVelocity * (
        upwindWeight * upwindConcentration + (1.0 - upwindWeight) * downwindConcentration
    )


##########################
# Define solution vector #
##########################

solution = np.zeros(gridGeometry.numDofs)


###################
# Enable plotting #
###################


@gridFunction(gridView)
def solutionGridFunction(element, x):
    elementIdx = elementMapper.index(element)
    return solution[elementIdx]


def plot(time):
    if plotting and dimension == 2:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        solutionGridFunction.plot(figure=fig, clim=[0, 1 + 1e-6], gridLines=None)
        ax.set_title("t = " + "{:0.2f}".format(time))
        plt.show()
        plt.pause(1e-3)
    gridView.writeVTK(
        problem.name + "-solution-{:0.2f}".format(time).replace(".", ""),
        celldata={"solution": solutionGridFunction},
        outputType=OutputType.ascii,
    )


#######################
# Initialize solution #
#######################

for element in gridView.elements:
    elementIdx = elementMapper.index(element)
    solution[elementIdx] = problem.initial(element)
plot(time=0)


####################
# Implement update #
####################


def assembleUpdate():
    update = np.zeros(gridGeometry.numDofs)

    for element in gridView.elements:
        fvGeometry = gridGeometry.localView
        fvGeometry.bind(element)

        for scvf in fvGeometry.scvfs:
            insideScvIdx = scvf.insideScvIdx
            insideConcentration = solution[insideScvIdx]

            if scvf.boundary:
                bndType = problem.boundaryTypes(element, scvf)
                if bndType.isDirichlet:
                    outsideConcentration = problem.dirichlet(element, scvf)[0]
                else:
                    raise Exception("Only Dirichlet BCs are implemented!")

            else:
                outsideScvIdx = scvf.outsideScvIdx
                outsideConcentration = solution[outsideScvIdx]

            flux = (
                advectiveFlux(insideConcentration, outsideConcentration, scvf.unitOuterNormal)
                * scvf.area
            )

            update[insideScvIdx] -= flux

        for scv in fvGeometry.scvs:
            update[scv.dofIndex] /= scv.volume

    return update


###################
# Start time loop #
###################

dt = 0.5 / (cells * velocity.two_norm)  # CFL condition
timeLoop = TimeLoop(startTime=0.0, dt=dt, endTime=1.0, verbose=True)
timeLoop.setPeriodicCheckPoint(0.25)

timeLoop.start()
while not timeLoop.finished:
    update = assembleUpdate()
    solution += timeLoop.timeStepSize * update
    timeLoop.advanceTimeStep()
    timeLoop.reportTimeStep()

    if timeLoop.isCheckPoint:
        plot(time=timeLoop.time)

timeLoop.finalize()
