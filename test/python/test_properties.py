#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""A model must not write into the list of parents it was given.

The discretization tag is deduced from the grid geometry and added to the parents, and the
names among them are resolved to TypeTag objects. Doing either in the caller's own list makes
one model's construction visible to the next: a second model built from the same list inherits
the first model's discretization as well as its own.
"""

from dune.grid import structuredGrid
from dumux.common import Model
from dumux.discretization import GridGeometry

gridView = structuredGrid([0, 0], [1, 1], [3, 3])
boxGeometry = GridGeometry(gridView, discMethod="box")
tpfaGeometry = GridGeometry(gridView, discMethod="cctpfa")

parents = ["OneP"]

boxModel = Model(inheritsFrom=parents, gridGeometry=boxGeometry)
if parents != ["OneP"]:
    raise Exception(f"building a model changed the list of parents it was given: {parents}")

tpfaModel = Model(inheritsFrom=parents, gridGeometry=tpfaGeometry)
if parents != ["OneP"]:
    raise Exception(f"building a second model changed the list of parents: {parents}")


def parentNames(model):
    return sorted(parent.name for parent in model.inheritsFrom)


if parentNames(boxModel) != ["BoxModel", "OneP"]:
    raise Exception(f"the box model inherits from {parentNames(boxModel)}")
if parentNames(tpfaModel) != ["CCTpfaModel", "OneP"]:
    raise Exception(f"the cell-centred model inherits from {parentNames(tpfaModel)}")

print("two models built from one list of parents keep their own discretizations")

description = boxModel.describe()
if description["typeTag"] != boxModel.name:
    raise Exception(f"describe() names the model {description['typeTag']}")
if sorted(description["inheritsFrom"]) != ["BoxModel", "OneP"]:
    raise Exception(f"describe() reports the parents as {description['inheritsFrom']}")
if description["properties"].get("Scalar") != "double":
    raise Exception(f"describe() reports the properties as {description['properties']}")
if description["gridGeometry"] != boxGeometry._typeName:
    raise Exception("describe() reports a grid geometry the model was not built on")

print("a model describes itself without generating its C++ header")
