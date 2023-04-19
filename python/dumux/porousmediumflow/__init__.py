# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Classes and functions related to the porousmedium flow models"""

import numpy as np
from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dumux.common import Property
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


@cppWrapperCreator
def _createPorousMediumFlowProblemDecorator(
    gridGeometry, spatialParams, enableInternalDirichletConstraints=False
):
    """A problem decorator generator for Python problems

    Usage:
        from dumux.common import PorousMediumFlowProblem
        @PorousMediumFlowProblem(gridGeometry)
        class MyProblem:
            ...
    """

    def createModule(numEq):
        priVars = f"Dune::FieldVector<double, {numEq}>"
        ggType = gridGeometry._typeName
        spType = spatialParams._typeName
        enableIDC = "true" if enableInternalDirichletConstraints else "false"
        problemType = (
            "Dumux::Python::PorousMediumFlowProblem" f"<{ggType}, {spType}, {priVars}, {enableIDC}>"
        )
        includes = [
            *(gridGeometry._includes),
            *(spatialParams._includes),
            *["dumux/python/porousmediumflow/problem.hh"],
        ]
        moduleName = "fvproblem_" + hashIt(problemType)
        generator = SimpleGenerator("PorousMediumFlowProblem", "Dumux::Python")
        module = generator.load(includes, problemType, moduleName, holder="std::shared_ptr")
        return module

    def decoratePorousMediumFlowProblem(cls):
        module = createModule(cls.numEq)

        def createPorousMediumFlowProblem():
            return module.PorousMediumFlowProblem(gridGeometry, spatialParams, cls())

        return createPorousMediumFlowProblem

    return decoratePorousMediumFlowProblem


@cppWrapperClassAlias(creator=_createPorousMediumFlowProblemDecorator)
class PorousMediumFlowProblem:
    """A class alias used to create a problem decorator Python problems"""


@cppWrapperCreator
def _createPorousMediumFlowVelocityOutput(*, gridVariables):
    """Create a PorousMediumFlowVelocityOutput"""

    includes = gridVariables._includes
    includes += ["dumux/python/porousmediumflow/velocityoutput.hh", "dumux/io/velocityoutput.hh"]
    fluxVarsType = (
        f"Dumux::GetPropType<{gridVariables.model.cppType}, Dumux::Properties::FluxVariables>"
    )
    typeName = f"Dumux::PorousMediumFlowVelocityOutput<{gridVariables._typeName}, {fluxVarsType}>"
    moduleName = "porousmediumflowvelocityoutput_" + hashIt(typeName)
    baseClass = [f"Dumux::VelocityOutput<{gridVariables._typeName}>"]
    generator = SimpleGenerator("PorousMediumFlowVelocityOutput", "Dumux::Python")
    module = generator.load(
        includes,
        typeName,
        moduleName,
        holder="std::shared_ptr",
        preamble=gridVariables.model.cppHeader,
        baseClasses=baseClass,
    )
    return module.PorousMediumFlowVelocityOutput(gridVariables)


@cppWrapperClassAlias(creator=_createPorousMediumFlowVelocityOutput)
class PorousMediumFlowVelocityOutput:
    """A class alias used to create PorousMediumFlowVelocityOutput instances"""


@cppWrapperCreator
def _createFVSpatialParamsOnePDecorator(gridGeometry):
    """A spatial params decorator generator for Python spatial params

    from dumux.common import FVSpatialParamsOneP
    FVSpatialParamsOneP(gridGeometry)
    class MySpatialParams:
        ...
    """
    dim = gridGeometry.gridView.dimensionworld
    scalar = Property.fromCppType("double")
    scalarType = scalar.cppType
    permType = f"Dune::FieldMatrix<{scalarType}, {dim}, {dim}>"

    def decorateFVSpatialParamsOneP(cls):
        ggType = gridGeometry._typeName
        spatialParamsType = f"Dumux::Python::FVSpatialParamsOneP<{ggType}, {permType}>"
        includes = gridGeometry._includes + ["dumux/python/porousmediumflow/spatialparams.hh"]
        includes += ["dumux/python/common/fvspatialparams.hh"]
        moduleName = "fvspatialparamsonep_" + hashIt(spatialParamsType)
        holderType = f"std::shared_ptr<{spatialParamsType}>"
        generator = SimpleGenerator("FVSpatialParamsOneP", "Dumux::Python")
        module = generator.load(includes, spatialParamsType, moduleName, options=[holderType])

        def maybeConvertScalarToMatrix(permeabilityValue):
            if isinstance(permeabilityValue, float):
                matrix = np.zeros(shape=(dim, dim))
                np.fill_diagonal(matrix, permeabilityValue)
                return matrix.tolist()

            return permeabilityValue

        class Permeability:
            """Permeability decorator to make sure permeability has correct type"""

            def __init__(self, permeabilityFunction):
                self.permeabilityFunction = permeabilityFunction

            def __call__(self, element, scv, elemSol):
                result = self.permeabilityFunction(element, scv, elemSol)
                return maybeConvertScalarToMatrix(result)

        class PermeabilityAtPos:
            """PermeabilityAtPos decorator to make sure permeability has correct type"""

            def __init__(self, permeabilityFunction):
                self.permeabilityFunction = permeabilityFunction

            def __call__(self, globalPos):
                result = self.permeabilityFunction(globalPos)
                return maybeConvertScalarToMatrix(result)

        def createFVSpatialParamsOneP():
            spatialParams = cls()
            if hasattr(cls, "permeability"):
                cls.permeability = Permeability(spatialParams.permeability)
            if hasattr(cls, "permeabilityAtPos"):
                cls.permeabilityAtPos = PermeabilityAtPos(spatialParams.permeabilityAtPos)
            return module.FVSpatialParamsOneP(gridGeometry, cls())

        return createFVSpatialParamsOneP

    return decorateFVSpatialParamsOneP


@cppWrapperClassAlias(creator=_createFVSpatialParamsOnePDecorator)
class FVSpatialParamsOneP:
    """Class alias used to decorate a Python Finite Volume Spatial Params"""
