from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

# A problem decorator generator for Python problems
#
# from dumux.common import PorousMediumFlowProblem
# @PorousMediumFlowProblem(gridGeometry)
# class MyProblem:
#    ...
#
def PorousMediumFlowProblem(gridGeometry, spatialParams, enableInternalDirichletConstraints=False):
    def createModule(numEq):
        priVarType = f"Dune::FieldVector<double, {numEq}>"
        ggType = gridGeometry._typeName
        spatialParamsType = spatialParams._typeName
        enableIntDirConstraint = "true" if enableInternalDirichletConstraints else "false"
        problemType = f"Dumux::Python::PorousMediumFlowProblem<{ggType}, {priVarType}, {spatialParamsType}, {enableIntDirConstraint}>"
        includes = (
            gridGeometry._includes
            + spatialParams._includes
            + ["dumux/python/porousmediumflow/problem.hh"]
        )
        moduleName = "fvproblem_" + hashIt(problemType)
        generator = SimpleGenerator("PorousMediumFlowProblem", "Dumux::Python")
        module = generator.load(includes, problemType, moduleName, holder="std::shared_ptr")
        return module

    def PorousMediumFlowProblemDecorator(Cls):
        module = createModule(Cls.numEq)

        def createPorousMediumFlowProblem():
            return module.PorousMediumFlowProblem(gridGeometry, spatialParams, Cls())

        return createPorousMediumFlowProblem

    return PorousMediumFlowProblemDecorator


def PorousMediumFlowVelocityOutput(*, gridVariables):
    includes = gridVariables._includes
    includes += ["dumux/python/porousmediumflow/velocityoutput.hh", "dumux/io/velocityoutput.hh"]
    fluxVarsType = (
        f"Dumux::GetPropType<{gridVariables._model.getTypeTag()}, Dumux::Properties::FluxVariables>"
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
        preamble=gridVariables._model.getProperties(),
        baseClasses=baseClass,
    )
    return module.PorousMediumFlowVelocityOutput(gridVariables)
