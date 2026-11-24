#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""The spatial-params decorator has to leave the decorated class alone.

The permeability callback is wrapped so that a scalar reaches C++ as a tensor. Wrapping it on
the class rather than on the instance would make the wrapper global: the object C++ holds and
the object whose callback it ends up calling would be different ones, so state read in the
callback would not be the state of the parameters that were handed over, and creating a second
set of parameters would change the answers the first one gives.
"""

from dune.grid import structuredGrid
from dumux.discretization import GridGeometry
from dumux.porousmediumflow import FVSpatialParamsOneP

gridView = structuredGrid([0, 0], [1, 1], [3, 3])
gridGeometry = GridGeometry(gridView, discMethod="box")

permeability = 1e-12


class SpatialParams:
    """Stateful on purpose: it reads what its permeability is when it is constructed."""

    instances = 0

    def __init__(self):
        SpatialParams.instances += 1
        self.permeabilityValue = permeability

    def permeability(self, element, scv, elemSol):
        return self.permeabilityValue

    def porosityAtPos(self, globalPos):
        return 0.4


originalCallback = SpatialParams.permeability
makeSpatialParams = FVSpatialParamsOneP(gridGeometry=gridGeometry)(SpatialParams)

first = makeSpatialParams()
if SpatialParams.permeability is not originalCallback:
    raise Exception("creating spatial params replaced the permeability callback on the class")
if SpatialParams.instances != 1:
    raise Exception(f"one set of spatial params constructed {SpatialParams.instances} instances, "
                    f"so the instance C++ holds is not the one the callback reads")

permeability = 1e-10
second = makeSpatialParams()
if SpatialParams.permeability is not originalCallback:
    raise Exception("creating a second set of spatial params replaced the class callback, which "
                    "would make the first set answer with the second one's state")
if SpatialParams.instances != 2:
    raise Exception(f"two sets of spatial params constructed {SpatialParams.instances} instances")

print("the spatial-params decorator wraps the instance, not the class")
