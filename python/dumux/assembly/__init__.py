# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# pylint: skip-file
# until that decorator bit can be removed again

"""Classes and function related to the assembly of linear systems"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


def decoratePre(pre):
    def wrappedPre(*args, **kwargs):
        preamble = pre(*args, **kwargs)
        newPreamble = ""
        for line in preamble.split("\n"):
            newPreamble += line + "\n"
            if line.startswith("#include <config.h>"):
                newPreamble += "#undef DUMUX_MULTITHREADING_BACKEND\n"
                newPreamble += "#define DUMUX_MULTITHREADING_BACKEND Serial\n"
        return newPreamble

    return wrappedPre


myAttributes = vars(SimpleGenerator).copy()
myAttributes["pre"] = decoratePre(myAttributes["pre"])
MySimpleGenerator = type("MySimpleGenerator", (object,), myAttributes)


@cppWrapperCreator
def _createFVAssembler(
    *, problem, gridVariables, model, diffMethod="numeric", isImplicit=True, multithreading=False
):
    """
    Create an FVAssembler object

    Args:
        problem: A problem instance (boundary & initial conditions)
        gridVariables: A grid variable instance (primary & secondary variables defined on a grid)
        model: A DuMux model configuration instance
        diffMethod (str): The method to compute derivatives of the residual (numeric or analytic)
        isImplicit (bool): If the time discretization method is implicit or explicit
        multithreading (bool): Keep the configured backend instead of forcing Serial. Only for a
            problem that enters no Python code during assembly; the default cannot.

    Returns:
        An assembler object (Python-bindings of the corresponding C++ type)

    Usage:
        assembler = FVAssembler(
            problem=problem, gridVariables=gridVars, model=model, diffMethod=diffMethod
        )
    """

    assemblerType = (
        f"Dumux::FVAssembler<{model.cppType}, {_cppDiffMethod(diffMethod)}, {int(isImplicit)}>"
    )
    return _loadAssembler(
        assemblerType,
        ["dumux/assembly/fvassembler.hh", "dumux/python/assembly/fvassembler.hh"],
        problem=problem,
        gridVariables=gridVariables,
        model=model,
        multithreading=multithreading,
    )


@cppWrapperClassAlias(creator=_createFVAssembler)
class FVAssembler:
    """Class alias used to instantiate an FVAssembler"""


@cppWrapperCreator
def _createAssembler(
    *, problem, gridVariables, model, diffMethod="numeric", isImplicit=True, multithreading=False
):
    """
    Create the assembler the model's type tag needs

    A model on finite-volume grid variables gets an FVAssembler, a model on the general grid
    variables the general assembler, which is the only one that can assemble it. The arguments are
    those of FVAssembler.
    """

    diff = _cppDiffMethod(diffMethod)
    assemblerType = f"Dumux::Python::AssemblerFor<{model.cppType}, {diff}, {int(isImplicit)}>"
    return _loadAssembler(
        assemblerType,
        ["dumux/python/assembly/assembler.hh"],
        problem=problem,
        gridVariables=gridVariables,
        model=model,
        multithreading=multithreading,
    )


@cppWrapperClassAlias(creator=_createAssembler)
class Assembler:
    """Class alias used to instantiate the assembler a model needs"""


def _cppDiffMethod(diffMethod):
    if diffMethod == "numeric":
        return "Dumux::DiffMethod::numeric"
    if diffMethod == "analytic":
        return "Dumux::DiffMethod::analytic"
    raise ValueError(f"Unknown diffMethod {diffMethod}")


def _loadAssembler(assemblerType, headers, *, problem, gridVariables, model, multithreading):
    includes = problem._includes + problem.gridGeometry()._includes + headers

    moduleName = "fvassembler_" + hashIt(assemblerType + str(multithreading))
    # remark: use SimpleGenerator again starting with dune 2.9
    generator = (SimpleGenerator if multithreading else MySimpleGenerator)(
        "FVAssembler", "Dumux::Python"
    )
    module = generator.load(
        includes,
        assemblerType,
        moduleName,
        holder="std::shared_ptr",
        preamble=model.cppHeader,
        # by default the assembler is compiled with the Serial backend, as a problem that calls
        # into Python during assembly is not thread-safe
        # the following is nicer but only works with dune > 2.8
        # extraCMake=[
        #     "target_compile_definitions(TARGET PUBLIC DUMUX_MULTITHREADING_BACKEND=Serial)"
        # ],
    )
    return module.FVAssembler(problem, problem.gridGeometry(), gridVariables)
