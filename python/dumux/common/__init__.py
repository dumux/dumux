"""
The DuMux common module
containing classes and functions needed for most simulations
"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

from dumux.common.properties import Model, Property
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias

from ._common import *


@cppWrapperCreator
def _createFVProblemDecorator(gridGeometry, enableInternalDirichletConstraints=False):
    """A problem decorator generator for Python problems

    from dumux.common import FVProblem
    @FVProblem(gridGeometry)
    class MyProblem:
        ...
    """

    def createModule(numEq):
        priVarType = f"Dune::FieldVector<double, {numEq}>"
        ggType = gridGeometry._typeName
        enableIntDirConstraint = "true" if enableInternalDirichletConstraints else "false"
        problemType = f"Dumux::Python::FVProblem<{ggType}, {priVarType}, {enableIntDirConstraint}>"
        includes = gridGeometry._includes + ["dumux/python/common/fvproblem.hh"]
        moduleName = "fvproblem_" + hashIt(problemType)
        holderType = f"std::shared_ptr<{problemType}>"
        generator = SimpleGenerator("FVProblem", "Dumux::Python")
        module = generator.load(includes, problemType, moduleName, options=[holderType])
        return module

    def decorateFVProblem(cls):
        module = createModule(cls.numEq)

        def createFVProblem():
            return module.FVProblem(gridGeometry, cls())

        return createFVProblem

    return decorateFVProblem


@cppWrapperClassAlias(creator=_createFVProblemDecorator)
class FVProblem:
    """Class alias used to decorate a Python finite volume problem"""


@cppWrapperCreator
def _createBoundaryTypes(numEq=1):
    """Create BoundaryTypes instances"""

    # only compile this once per numEq
    cacheKey = f"BoundaryTypes_{numEq}"
    try:
        return globals()[cacheKey]()
    except KeyError:
        includes = ["dumux/python/common/boundarytypes.hh"]
        typeName = f"Dumux::BoundaryTypes<{numEq}>"
        moduleName = "boundarytypes_" + hashIt(typeName)
        generator = SimpleGenerator("BoundaryTypes", "Dumux::Python")
        module = generator.load(includes, typeName, moduleName)
        globals().update({cacheKey: module.BoundaryTypes})
    return globals()[cacheKey]()


@cppWrapperClassAlias(creator=_createBoundaryTypes)
class BoundaryTypes:
    """Class alias used to create a BoundaryTypes instance"""
