#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Grid variables and assembler of a model on the general grid variables

The model is the hyperelastic one, which declares the general grid variables, so `GridVariables`
and `Assembler` have to pick the general grid variables and the general assembler. The problem
imposes a homogeneous deformation on the whole boundary, whose discrete solution is exact on a box
grid: Newton has to reach u = (F - I)x to round-off.
"""

import sys

import numpy as np

from dune.common.hashit import hashIt
from dune.generator.generator import SimpleGenerator
from dune.grid import structuredGrid
from dune.common import FieldVector
from dune.istl import BCRSMatrix, BlockVector

from dumux.common import initParameters, Model, Property
from dumux.common.properties import CppTypeTag
from dumux.discretization import GridGeometry, GridVariables
from dumux.assembly import Assembler

PROBLEM_HEADER = "test/python/test_generalgridvariables.hh"

initParameters(
    argv=sys.argv,
    params={"SpatialParams.YoungsModulus": 1.0, "SpatialParams.PoissonRatio": 0.3},
)


def homogeneousDeformationProblem(model, gridGeometry, F):
    typeName = f"Dumux::HomogeneousDeformationProblem<{model.cppType}>"
    generator = SimpleGenerator("HomogeneousDeformationProblem", "Dumux::Python")
    module = generator.load(
        gridGeometry._includes + [PROBLEM_HEADER],
        typeName,
        "homogeneousdeformation_" + hashIt(typeName),
        holder="std::shared_ptr",
        preamble=model.cppHeader,
    )
    return module.HomogeneousDeformationProblem(gridGeometry, F)


# a stretch, a shear and a rotation at once
F = [[1.1, 0.3], [-0.2, 0.9]]

gridView = structuredGrid([0.0, 0.0], [1.0, 1.0], [6, 6])
gridGeometry = GridGeometry(gridView=gridView, discMethod="box")

hyperelastic = CppTypeTag(name="Hyperelastic", include="dumux/solidmechanics/hyperelastic/model.hh")
model = Model(inheritsFrom=[hyperelastic], gridGeometry=gridGeometry)
model["Problem"] = Property.fromCppType(
    "HomogeneousDeformationProblem<TypeTag>", cppIncludes=gridGeometry._includes + [PROBLEM_HEADER]
)
# cached, so the assembler reads the variables `update` writes, not ones computed from the solution
model["EnableGridVolumeVariablesCache"] = Property.fromValue(True)
problem = homogeneousDeformationProblem(model, gridGeometry, F)

gridVariables = GridVariables(problem=problem, model=model)
if gridVariables.gridDiscretization.numDofs != gridGeometry.numDofs:
    raise Exception("The general grid variables do not hold the problem's grid discretization")


def displacementVector():
    return BlockVector(FieldVector([0.0, 0.0]), gridGeometry.numDofs)


dofPosition = {}
for element in gridView.elements:
    fvGeometry = gridGeometry.boundLocalView(element)
    for scv in fvGeometry.scvs:
        dofPosition[scv.dofIndex] = scv.dofPosition

# a start that is not a homogeneous deformation, so the tangent varies in space and Newton has to
# iterate rather than land on the solution in one step
sol = displacementVector()
for dof, x in dofPosition.items():
    bump = 0.1 * np.sin(np.pi * x[0]) * np.sin(np.pi * x[1])
    sol[dof] = [bump, -bump]
gridVariables.init(sol)
assembler = Assembler(
    problem=problem, gridVariables=gridVariables, model=model, diffMethod="numeric"
)


# the assembler hands its Jacobian to Python, which needs the 2x2-block matrix type registered
BCRSMatrix((2, 2))


def toArray(vector):
    return np.array([[vector[i][0], vector[i][1]] for i in range(len(vector))])


def toDense(matrix):
    dense = np.zeros((2 * matrix.rows, 2 * matrix.rows))
    for rowIndex, row in enumerate(matrix):
        for columnIndex, block in row.enumerate:
            dense[2 * rowIndex : 2 * rowIndex + 2, 2 * columnIndex : 2 * columnIndex + 2] = (
                np.asarray(block, float).reshape(2, 2)
            )
    return dense


residualNorms = []
for _ in range(10):
    assembler.assembleJacobianAndResidual(sol)
    residual = toArray(assembler.residual)
    residualNorms.append(np.max(np.abs(residual)))
    if residualNorms[-1] < 1e-12:
        break

    update = np.linalg.solve(toDense(assembler.jacobian), residual.ravel()).reshape(-1, 2)
    for i, (du0, du1) in enumerate(update):
        sol[i] = [sol[i][0] - du0, sol[i][1] - du1]
    gridVariables.update(sol)

print("residual max-norm per Newton iteration:", " ".join(f"{r:.2e}" for r in residualNorms))
if residualNorms[-1] >= 1e-12:
    raise Exception("Newton did not converge")

if len(residualNorms) < 4:
    raise Exception("Newton converged too fast for the start to have tested repeated assembly")

error = 0.0
for dof, x in dofPosition.items():
    for i in range(2):
        exact = sum(F[i][j] * x[j] for j in range(2)) - x[i]
        error = max(error, abs(sol[dof][i] - exact))

print(f"max displacement error: {error:.2e}")
if error > 1e-10:
    raise Exception(f"The homogeneous deformation is not reproduced, error {error:.2e}")
