# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Boundary types generator
"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


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
