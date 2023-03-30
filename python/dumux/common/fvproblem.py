# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Finite volume problem generator
"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


@cppWrapperCreator
def _createFVProblemDecorator(
    gridGeometry, spatialParams=None, enableInternalDirichletConstraints=False
):
    """A problem decorator generator for Python problems

    from dumux.common import FVProblem
    @FVProblem(gridGeometry)
    class MyProblem:
        ...
    """

    def createModule(numEq):
        ggType = gridGeometry._typeName
        if spatialParams is not None:
            spType = spatialParams._typeName
        else:
            spType = f"Dumux::Python::FVSpatialParams<{ggType}>"

        priVarType = f"Dune::FieldVector<double, {numEq}>"
        enableIntDirConstraint = "true" if enableInternalDirichletConstraints else "false"
        problemType = (
            f"Dumux::Python::FVProblem<{ggType}, {spType}, {priVarType}, {enableIntDirConstraint}>"
        )
        includes = gridGeometry._includes + ["dumux/python/common/fvproblem.hh"]
        moduleName = "fvproblem_" + hashIt(problemType)
        holderType = f"std::shared_ptr<{problemType}>"
        generator = SimpleGenerator("FVProblem", "Dumux::Python")
        module = generator.load(includes, problemType, moduleName, options=[holderType])
        return module

    def decorateFVProblem(cls):
        module = createModule(cls.numEq)

        def createFVProblem():
            if spatialParams is not None:
                return module.FVProblem(gridGeometry, spatialParams, cls())

            return module.FVProblem(gridGeometry, cls())

        return createFVProblem

    return decorateFVProblem


@cppWrapperClassAlias(creator=_createFVProblemDecorator)
class FVProblem:
    """Class alias used to decorate a Python finite volume problem"""
