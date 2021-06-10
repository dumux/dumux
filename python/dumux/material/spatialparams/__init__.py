from dune.generator.generator import SimpleGenerator
from dumux.common import Property
import numpy as np


def OnePSpatialParams(*, gridGeometry, scalar=None):
    moduleName = "spatialparams"
    dim = gridGeometry.gridView.dimension
    includes = gridGeometry._includes + [
        "dumux/python/material/spatialparams/spatialparams.hh",
        "dune/common/fmatrix.hh",
    ]

    if scalar is None:
        scalar = Property(type="double")

    typeName = f"Dumux::Python::FVSpatialParamsOneP<{gridGeometry._typeName}, {scalar._typeName}, Dune::FieldMatrix<{scalar._typeName}, {dim}, {dim}>>"

    def OnePSpatialParamsDecorator(Cls):
        generator = SimpleGenerator("OnePSpatialParams", "Dumux::Python")
        module = generator.load(includes, typeName, moduleName, holder="std::shared_ptr")

        def maybeConvertScalarToMatrix(permeabilityValue):
            if isinstance(permeabilityValue, float):
                matrix = np.zeros(shape=(dim, dim))
                np.fill_diagonal(matrix, permeabilityValue)
                return matrix.tolist()
            else:
                return permeabilityValue

        class Permeability:
            def __init__(self, permeabilityFunction):
                self.permeabilityFunction = permeabilityFunction

            def __call__(self, element, scv, elemSol):
                result = self.permeabilityFunction(element, scv, elemSol)
                return maybeConvertScalarToMatrix(result)

        class PermeabilityAtPos:
            def __init__(self, permeabilityFunction):
                self.permeabilityFunction = permeabilityFunction

            def __call__(self, globalPos):
                result = self.permeabilityFunction(globalPos)
                return maybeConvertScalarToMatrix(result)

        def createSpatialParams():
            sp = Cls()
            if hasattr(Cls, "permeability"):
                Cls.permeability = Permeability(sp.permeability)
            else:
                Cls.permeabilityAtPos = PermeabilityAtPos(sp.permeabilityAtPos)
            spatialParams = module.OnePSpatialParams(gridGeometry, sp)
            return spatialParams

        return createSpatialParams

    return OnePSpatialParamsDecorator
