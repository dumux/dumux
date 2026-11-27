#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The instationary assembly through the FVAssembler bindings

A one-phase ideal gas on a closed domain: the instationary residual differs from the stationary
one by exactly the storage term, porosity * volume * (rho(p) - rho(pPrev)) / dt in each cell,
so it shows that the assembler uses the time step size and the previous solution it was given.
The assembler stores only the address of the previous solution, so the binding has to keep that
vector alive, and release it once it is replaced.
"""

import gc
import sys
import weakref

from dune.grid import structuredGrid
from dune.istl import blockVector

from dumux.common import initParameters, BoundaryTypes, Model, Property, TimeLoop
from dumux.discretization import GridGeometry, GridVariables
from dumux.assembly import FVAssembler
from dumux.porousmediumflow import PorousMediumFlowProblem, FVSpatialParamsOneP
from dumux.material import FluidSystem, Component

initParameters(argv=sys.argv, params={"Problem.EnableGravity": False})

POROSITY = 0.4
TEMPERATURE = 293.15
CELLS = 4
DT = 2.0

gridView = structuredGrid([0, 0], [1, 1], [CELLS, CELLS])
gridGeometry = GridGeometry(gridView=gridView, discMethod="cctpfa")

model = Model(inheritsFrom=["OneP"], gridGeometry=gridGeometry)
air = Component("Air")
model["FluidSystem"] = Property.fromInstance(
    FluidSystem("OnePGas", component=air, scalar=model["Scalar"])
)


@FVSpatialParamsOneP(gridGeometry=gridGeometry)
class SpatialParams:
    def permeabilityAtPos(self, globalPos):
        return 1e-12

    def porosityAtPos(self, globalPos):
        return POROSITY

    def temperatureAtPos(self, globalPos):
        return TEMPERATURE


spatialParams = SpatialParams()
model["SpatialParams"] = Property.fromInstance(spatialParams)


@PorousMediumFlowProblem(gridGeometry=gridGeometry, spatialParams=spatialParams)
class Problem:
    numEq = 1

    def boundaryTypesAtPos(self, globalPos):
        bTypes = BoundaryTypes(self.numEq)
        bTypes.setNeumann()
        return bTypes

    def neumann(self, element, fvGeometry, scvf):
        return 0.0

    def sourceAtPos(self, globalPos):
        return 0.0

    def name(self):
        return "fvassembler"

    def paramGroup(self):
        return ""


problem = Problem()
model["Problem"] = Property.fromInstance(problem)

numDofs = gridGeometry.numDofs
pressure = [1.0e5 + 2.0e3 * i for i in range(numDofs)]
previousPressure = [p - 5.0e2 * (1 + i % 3) for i, p in enumerate(pressure)]


def vector(values):
    result = blockVector(numDofs)
    for i, value in enumerate(values):
        result[i] = value
    return result


sol = vector(pressure)
gridVars = GridVariables(problem=problem, model=model)
gridVars.init(sol)
assembler = FVAssembler(problem=problem, gridVariables=gridVars, model=model)

if not assembler.isStationaryProblem:
    raise Exception("An assembler without a time loop has to assemble the stationary problem")
assembler.assembleJacobianAndResidual(sol)
stationaryResidual = [r[0] for r in assembler.residual]

assembler.setTimeLoop(TimeLoop(0.0, DT, DT, False))
if assembler.isStationaryProblem:
    raise Exception("An assembler with a time loop has to assemble the instationary problem")

previous = vector(previousPressure)
previousAlive = weakref.ref(previous)
assembler.setPreviousSolution(previous)
del previous
gc.collect()
if previousAlive() is None:
    raise Exception("The previous solution was freed while the assembler still refers to it")

assembler.assembleJacobianAndResidual(sol)
cellVolume = 1.0 / CELLS**2
maxError = 0.0
for i, r in enumerate(assembler.residual):
    density = air.gasDensity(TEMPERATURE, pressure[i])
    previousDensity = air.gasDensity(TEMPERATURE, previousPressure[i])
    storage = POROSITY * cellVolume * (density - previousDensity) / DT
    maxError = max(maxError, abs((r[0] - stationaryResidual[i]) - storage) / abs(storage))
print(f"largest relative error of the storage term: {maxError:.2e}")
if maxError > 1e-10:
    raise Exception("The instationary residual is not the stationary one plus the storage term")

assembler.setPreviousSolution(vector(pressure))
gc.collect()
if previousAlive() is not None:
    raise Exception("A replaced previous solution is still held by the assembler")
