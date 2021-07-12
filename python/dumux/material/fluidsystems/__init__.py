from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

# construct a FluidSystem
# the FluidSystem is JIT compiled
def FluidSystem(*, type, component, scalar="double"):
    print(component._includes)
    includes = component._includes + ["dumux/python/material/fluidsystems/fluidsystem.hh"]

    availabeFluidSystems = {"OnePLiquid" : [["dumux/material/fluidsystems/1pliquid.hh"],
                                             "Dumux::FluidSystems::OnePLiquid<{}, {}>".format(scalar._typeName, component._typeName)],
                            "OnePGas" :    [["dumux/material/fluidsystems/1pgas.hh"],
                                             "Dumux::FluidSystems::OnePGas<{}, {}>".format(scalar._typeName, component._typeName)]}

    if type not in availabeFluidSystems:
        raise NotImplementedError("FluidSystem of type " + type + " not implemented.\n"
                                  "Availabe types are " + ", ".join(availabeFluidSystems.keys()))

    includes += availabeFluidSystems[type][0]
    typeName = availabeFluidSystems[type][1]

    moduleName = "FluidSystem_" + hashIt(typeName)
    generator = SimpleGenerator("FluidSystem", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName)
    return module.FluidSystem()
