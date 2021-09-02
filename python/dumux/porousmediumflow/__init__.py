"""Classes and functions related to the porousmedium flow models"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
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
            "Dumux::Python::PorousMediumFlowProblem" f"<{ggType}, {priVars}, {spType}, {enableIDC}>"
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
