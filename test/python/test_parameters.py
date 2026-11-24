#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Runtime parameters have to survive the round trip through the parameter tree.

The tree holds text, so a Python value is written out and read back. Two spellings differ
between the languages and have to be translated in both directions: a truth value (``True`` is
not ``true``, and every non-empty string, ``"false"`` included, is true in Python) and a list
(``[10, 12]`` is not ``10 12``).
"""

import sys

from dumux.common import getParam, initParameters

initParameters(
    argv=sys.argv[:1],
    params={
        "Problem.EnableGravity": False,
        "Problem.Name": "roundtrip",
        "Grid.Cells": [10, 12],
        "Grid.UpperRight": [1.5, 2.5],
        "Newton.MaxSteps": 8,
        "TimeLoop.DtInitial": 0.5,
    },
)


def check(name, got, expected):
    if type(got) is not type(expected) or got != expected:
        raise Exception(f"{name} read back as {got!r}, expected {expected!r}")
    print(f"  {name}: {got!r}")


check("Problem.EnableGravity", getParam("Problem.EnableGravity", default=True), False)
check("Problem.Name", getParam("Problem.Name"), "roundtrip")
check("Problem.Name", getParam("Problem.Name", default="other"), "roundtrip")
check("Grid.Cells", getParam("Grid.Cells", default=[1, 1]), [10, 12])
check("Grid.UpperRight", getParam("Grid.UpperRight", default=[1.0, 1.0]), [1.5, 2.5])
check("Newton.MaxSteps", getParam("Newton.MaxSteps", default=1), 8)
check("TimeLoop.DtInitial", getParam("TimeLoop.DtInitial", default=1.0), 0.5)

check("Unset.Flag", getParam("Unset.Flag", default=True), True)
check("Unset.Cells", getParam("Unset.Cells", default=[3, 4]), [3, 4])
check("Unset.Name", getParam("Unset.Name", default="fallback"), "fallback")

print("parameters keep their type across the parameter tree")
