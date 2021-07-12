from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

def FVAssembler(*, problem, gridVariables, model, diffMethod='numeric', isImplicit=True):

    TypeTag = model.getTypeTag()

    if diffMethod == 'numeric':
        dm = 'Dumux::DiffMethod::numeric'
    elif diffMethod == 'analytic':
        dm = 'Dumux::DiffMethod::analytic'
    else:
        raise ValueError("Unknown diffMethod {}".format(diffMethod))

    assemblerType = "Dumux::FVAssembler<{}, {}>".format(TypeTag, dm, int(isImplicit))
    includes = problem._includes + problem.gridGeometry()._includes + ["dumux/assembly/fvassembler.hh"]
    includes += ["dumux/python/assembly/fvassembler.hh"]

    moduleName = "fvassembler_" + hashIt(assemblerType)
    holderType = "std::shared_ptr<{}>".format(assemblerType)
    generator = SimpleGenerator("FVAssembler", "Dumux::Python")
    module = generator.load(includes, assemblerType, moduleName, options=[holderType], preamble=model.getProperties())
    return module.FVAssembler(problem, problem.gridGeometry(), gridVariables)
