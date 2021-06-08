from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt


def FVAssembler(*, problem, gridVariables, model, diffMethod="numeric", isImplicit=True):

    TypeTag = model.getTypeTag()

    if diffMethod == "numeric":
        dm = "Dumux::DiffMethod::numeric"
    elif diffMethod == "analytic":
        dm = "Dumux::DiffMethod::analytic"
    else:
        raise ValueError(f"Unknown diffMethod {diffMethod}")

    assemblerType = f"Dumux::FVAssembler<{TypeTag}, {dm}, {int(isImplicit)}>"
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
        preamble=model.getProperties(),
    )
    return module.FVAssembler(problem, problem.gridGeometry(), gridVariables)
