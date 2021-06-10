from ._common import *

from dumux.common.properties import Model, Property

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt


# A problem decorator generator for Python problems
#
# from dumux.common import FVProblem
# @FVProblem(gridGeometry)
# class MyProblem:
#    ...
#
def FVProblem(gridGeometry, spatialParams):

    def createModule(numEq):
        priVarType = "Dune::FieldVector<double, {}>".format(numEq)
        ggType = gridGeometry._typeName
        spatialParamsType = spatialParams._typeName
        problemType = "Dumux::Python::FVProblem<{}, {}, {}>".format(ggType, priVarType, spatialParamsType)
        includes = gridGeometry._includes + spatialParams._includes + ["dumux/python/common/fvproblem.hh"]
        moduleName = "fvproblem_" + hashIt(problemType)
        holderType = "std::shared_ptr<{}>".format(problemType)
        generator = SimpleGenerator("FVProblem", "Dumux::Python")
        module = generator.load(includes, problemType, moduleName, options=[holderType])
        return module

    def FVProblemDecorator(Cls):
        module = createModule(Cls.numEq)

        def createFVProblem():
            return module.FVProblem(gridGeometry, spatialParams, Cls())
        return createFVProblem

    return FVProblemDecorator


# Function for JIT copmilation of Dumux::BoundaryTypes
def BoundaryTypes(numEq=1):
    # only compile this once per numEq
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


def Parameters(*, file=None, dict={}):
    parametersType = "Dumux::Parameters"
    includes = ["dumux/common/parameters.hh", "dumux/python/common/parameters.hh"]
    moduleName = "parameters_" + hashIt(parametersType)
    generator = SimpleGenerator("Parameters", "Dumux::Python")
    module = generator.load(includes, parametersType, moduleName)

    # make sure all dict keys are strings
    for key in dict:
        if not isinstance(dict[key], str):
            dict[key] = str(dict[key])

    if file is not None:
        return module.Parameters(file, dict)
    else:
        return module.Parameters(dict)
