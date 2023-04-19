# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Fluid systems"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dumux.common import Property
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


@cppWrapperCreator
def _createFluidSystem(name, *, component, scalar: Property = None):
    """Construct a FluidSystem"""

    includes = component._includes + ["dumux/python/material/fluidsystems/fluidsystem.hh"]

    if scalar is None:
        scalar = Property.fromCppType("double")
    scalarType = scalar.cppType

    availableFluidSystems = {
        "OnePLiquid": {
            "includes": ["dumux/material/fluidsystems/1pliquid.hh"],
            "type": f"Dumux::FluidSystems::OnePLiquid<{scalarType}, {component._typeName}>",
        },
        "OnePGas": {
            "includes": ["dumux/material/fluidsystems/1pgas.hh"],
            "type": f"Dumux::FluidSystems::OnePGas<{scalarType}, {component._typeName}>",
        },
    }
    if name not in availableFluidSystems:
        raise NotImplementedError(
            "FluidSystem of type " + name + " not implemented.\n"
            "Available types are " + ", ".join(availableFluidSystems.keys())
        )

    includes += availableFluidSystems[name]["includes"]
    typeName = availableFluidSystems[name]["type"]

    moduleName = "fluidsystem_" + hashIt(typeName)
    generator = SimpleGenerator("FluidSystem", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName)
    return module.FluidSystem()


@cppWrapperClassAlias(creator=_createFluidSystem)
class FluidSystem:
    """Class alias used to instantiate a FluidSystem"""
