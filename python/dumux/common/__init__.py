from ._common import *

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

# A problem decorator generator for Python problems
#
# from dumux.common import FVProblem
# @FVProblem(gridGeometry)
# class MyProblem:
#    ...
#
def FVProblem(gridGeometry):

    def createModule(numEq):
        priVarType = "Dune::FieldVector<double, {}>".format(numEq)
        ggType = gridGeometry._typeName
        problemType = "Dumux::Python::FVProblem<{}, {}>".format(ggType, priVarType)
        includes = gridGeometry._includes + ["dumux/python/common/fvproblem.hh"]
        moduleName = "fvproblem_" + hashIt(problemType)
        holderType = "std::shared_ptr<{}>".format(problemType)
        generator = SimpleGenerator("FVProblem", "Dumux::Python")
        module = generator.load(includes, problemType, moduleName, options=[holderType])
        return module

    def FVProblemDecorator(Cls):
        module = createModule(Cls.numEq)
        def createFVProblem():
            return module.FVProblem(gridGeometry, Cls())
        return createFVProblem

    return FVProblemDecorator


# Function for JIT copmilation of Dumux::BoundaryTypes
def BoundaryTypes(numEq=1):
    # only copmile this once per numEq
    cacheKey = "BoundaryTypes_{}".format(numEq)
    try:
        return globals()[cacheKey]()
    except KeyError:
        includes = ["dumux/python/common/boundarytypes.hh"]
        typeName = "Dumux::BoundaryTypes<{}>".format(numEq)
        moduleName = "boundarytypes_" + hashIt(typeName)
        generator = SimpleGenerator("BoundaryTypes", "Dumux::Python")
        module = generator.load(includes, typeName, moduleName)
        globals().update({cacheKey : module.BoundaryTypes})
    return globals()[cacheKey]()
