# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Finite volume spatial params generator
"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


@cppWrapperCreator
def _createFVSpatialParamsDecorator(gridGeometry):
    """A spatial params decorator generator for Python spatial params

    from dumux.common import FVSpatialParams
    @FVSpatialParams(gridGeometry)
    class MySpatialParams:
        ...
    """

    def createModule():
        ggType = gridGeometry._typeName
        spatialParamsType = f"Dumux::Python::FVSpatialParams<{ggType}>"
        includes = gridGeometry._includes + ["dumux/python/common/fvspatialparams.hh"]
        moduleName = "fvspatialparams_" + hashIt(spatialParamsType)
        holderType = f"std::shared_ptr<{spatialParamsType}>"
        generator = SimpleGenerator("FVSpatialParams", "Dumux::Python")
        module = generator.load(includes, spatialParamsType, moduleName, options=[holderType])
        return module

    def decorateFVSpatialParams(cls):
        module = createModule()

        def createFVSpatialParams():
            return module.FVSpatialParams(gridGeometry, cls())

        return createFVSpatialParams

    return decorateFVSpatialParams


@cppWrapperClassAlias(creator=_createFVSpatialParamsDecorator)
class FVSpatialParams:
    """Class alias used to decorate a Python Finite Volume Spatial Params"""
