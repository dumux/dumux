from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

# construct a VtkOutputModule
# the grid geometry is JIT compiled
def VtkOutputModule(*, gridVariables, solutionVector, name):
    includes = gridVariables._includes + solutionVector._includes
    includes += ["dumux/python/io/vtkoutputmodule.hh", "dumux/io/vtkoutputmodule.hh"]
    typeName = "Dumux::VtkOutputModule<{}, {}>".format(
        gridVariables._typeName, solutionVector._typeName
    )
    moduleName = "vtkoutputmodule_" + hashIt(typeName)
    generator = SimpleGenerator("VtkOutputModule", "Dumux::Python")
    module = generator.load(
        includes, typeName, moduleName, preamble=gridVariables._model.getProperties()
    )
    return module.VtkOutputModule(gridVariables, solutionVector, name)
