#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""A problem class is checked where it is written, not where it is called.

C++ looks the callbacks up by name while it assembles, so a name that is missing or a signature
that does not fit is otherwise reported from inside an element loop, after the generated code
has been compiled. These checks run at decoration time and name the class and the callback.
"""

from dune.grid import structuredGrid
from dumux.common import FVProblem
from dumux.common.validation import CallbackError
from dumux.discretization import GridGeometry

gridView = structuredGrid([0, 0], [1, 1], [3, 3])
gridGeometry = GridGeometry(gridView, discMethod="box")


def refused(description, cls):
    """Decorate a class that should not be accepted, and report why it was not"""

    try:
        FVProblem(gridGeometry=gridGeometry)(cls)
    except CallbackError as error:
        print(f"  {description}: {error}")
        return
    raise Exception(f"{description}: accepted a class that cannot be assembled")


class NoNumEq:
    def boundaryTypesAtPos(self, globalPos):
        return None


class NumEqIsNotACount:
    numEq = 0

    def boundaryTypesAtPos(self, globalPos):
        return None


class NoBoundaryTypes:
    numEq = 1

    def dirichletAtPos(self, globalPos):
        return 0.0


class BoundaryTypesTakesTheWrongArguments:
    numEq = 1

    def boundaryTypesAtPos(self, x, y):
        return None


class NeumannForgotTheGeometry:
    numEq = 1

    def boundaryTypesAtPos(self, globalPos):
        return None

    def neumann(self, element, scvf):
        return 0.0


refused("a problem without numEq", NoNumEq)
refused("a problem with numEq = 0", NumEqIsNotACount)
refused("a problem that says nothing about its boundary", NoBoundaryTypes)
refused("a position callback taking a position apart", BoundaryTypesTakesTheWrongArguments)
refused("a face callback missing the element geometry", NeumannForgotTheGeometry)

print("a problem that cannot be assembled is refused where it is written")
