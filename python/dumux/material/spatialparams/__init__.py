"""Spatial parameter classes"""

import numpy as np
from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dumux.common import Property
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


@cppWrapperCreator
def _createOnePSpatialParamsDecorator(*, gridGeometry, scalar: Property = None):
    """Turn a Python spatial parameter class into an C++/Python hybrid class"""

    dim = gridGeometry.gridView.dimensionworld
    includes = gridGeometry._includes + [
        "dumux/python/material/spatialparams/spatialparams.hh",
        "dune/common/fmatrix.hh",
    ]

    if scalar is None:
        scalar = Property.fromCppType("double")
    scalarType = scalar.cppType

    permType = f"Dune::FieldMatrix<{scalarType}, {dim}, {dim}>"
    typeName = (
        f"Dumux::Python::FVSpatialParamsOneP<{gridGeometry._typeName}, {scalarType}, {permType}>"
    )
    moduleName = f"spatialparams_{hashIt(typeName)}"

    def decorateOnePSpatialParams(cls):
        generator = SimpleGenerator("OnePSpatialParams", "Dumux::Python")
        module = generator.load(includes, typeName, moduleName, holder="std::shared_ptr")

        def maybeConvertScalarToMatrix(permeabilityValue):
            if isinstance(permeabilityValue, float):
                matrix = np.zeros(shape=(dim, dim))
                np.fill_diagonal(matrix, permeabilityValue)
                return matrix.tolist()

            return permeabilityValue

        class Permeability:
            """Permeability decorator to make sure permeability has correct type"""

            def __init__(self, permeabilityFunction):
                self.permeabilityFunction = permeabilityFunction

            def __call__(self, element, scv, elemSol):
                result = self.permeabilityFunction(element, scv, elemSol)
                return maybeConvertScalarToMatrix(result)

        class PermeabilityAtPos:
            """PermeabilityAtPos decorator to make sure permeability has correct type"""

            def __init__(self, permeabilityFunction):
                self.permeabilityFunction = permeabilityFunction

            def __call__(self, globalPos):
                result = self.permeabilityFunction(globalPos)
                return maybeConvertScalarToMatrix(result)

        def createSpatialParams():
            spatialParams = cls()
            if hasattr(cls, "permeability"):
                cls.permeability = Permeability(spatialParams.permeability)
            if hasattr(cls, "permeabilityAtPos"):
                cls.permeabilityAtPos = PermeabilityAtPos(spatialParams.permeabilityAtPos)
            return module.OnePSpatialParams(gridGeometry, spatialParams)

        return createSpatialParams

    return decorateOnePSpatialParams


@cppWrapperClassAlias(creator=_createOnePSpatialParamsDecorator)
class OnePSpatialParams:
    """Class alias used to decorate Python spatial parameter implementations"""
