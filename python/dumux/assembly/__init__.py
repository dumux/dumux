"""Classes and function related to the assembly of linear systems"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


@cppWrapperCreator
def _createFVAssembler(*, problem, gridVariables, model, diffMethod="numeric", isImplicit=True):
    """
    Create an FVAssembler object

    Args:
        problem: A problem instance (boundary & initial conditions)
        gridVariables: A grid variable instance (primary & secondary variables defined on a grid)
        model: A DuMux model configuration instance
        diffMethod (str): The method to compute derivatives of the residual (numeric or analytic)
        isImplicit (bool): If the time discretization method is implicit or explicit

    Returns:
        An assembler object (Python-bindings of the corresponding C++ type)

    Usage:
        assembler = FVAssembler(
            problem=problem, gridVariables=gridVars, model=model, diffMethod=diffMethod
        )
    """

    if diffMethod == "numeric":
        cppDiffMethod = "Dumux::DiffMethod::numeric"
    elif diffMethod == "analytic":
        cppDiffMethod = "Dumux::DiffMethod::analytic"
    else:
        raise ValueError(f"Unknown diffMethod {diffMethod}")

    assemblerType = f"Dumux::FVAssembler<{model.cppType}, {cppDiffMethod}, {int(isImplicit)}>"
    includes = (
        problem._includes + problem.gridGeometry()._includes + ["dumux/assembly/fvassembler.hh"]
    )
    includes += ["dumux/python/assembly/fvassembler.hh"]

    moduleName = "fvassembler_" + hashIt(assemblerType)
    generator = SimpleGenerator("FVAssembler", "Dumux::Python")
    module = generator.load(
        includes,
        assemblerType,
        moduleName,
        holder="std::shared_ptr",
        preamble=model.cppHeader,
    )
    return module.FVAssembler(problem, problem.gridGeometry(), gridVariables)


@cppWrapperClassAlias(creator=_createFVAssembler)
class FVAssembler:
    """Class alias used to instantiate an FVAssembler"""
