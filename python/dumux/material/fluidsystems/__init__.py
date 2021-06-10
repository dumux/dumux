from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

# construct a FluidSystem
# the FluidSystem is JIT compiled
def FluidSystem(*, type, component, scalar="double"):
    includes = component._includes + ["dumux/python/material/fluidsystems/fluidsystem.hh"]

    availableFluidSystems = {
        "OnePLiquid": {
            "includes": ["dumux/material/fluidsystems/1pliquid.hh"],
            "type": f"Dumux::FluidSystems::OnePLiquid<{scalar._typeName}, {component._typeName}>",
        },
        "OnePGas": {
            "includes": ["dumux/material/fluidsystems/1pgas.hh"],
            "type": f"Dumux::FluidSystems::OnePGas<{scalar._typeName}, {component._typeName}>",
        },
    }
    if type not in availableFluidSystems:
        raise NotImplementedError(
            "FluidSystem of type " + type + " not implemented.\n"
            "Available types are " + ", ".join(availableFluidSystems.keys())
        )

    includes += availableFluidSystems[type]["includes"]
    typeName = availableFluidSystems[type]["type"]

    moduleName = "FluidSystem_" + hashIt(typeName)
    generator = SimpleGenerator("FluidSystem", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName)
    return module.FluidSystem()
