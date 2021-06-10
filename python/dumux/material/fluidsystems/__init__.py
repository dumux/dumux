from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

# construct a GridGeometry from a gridView
# the grid geometry is JIT compiled
def FluidSystem(scalar, component):
    print(component._includes)
    includes = component._includes + ["dumux/python/material/fluidsystems/fluidsystem.hh"]
    includes += ["dumux/material/fluidsystems/1pliquid.hh"]


    typeName = "Dumux::FluidSystems::OnePLiquid<{}, {}>".format(scalar._typeName, component._typeName)


    moduleName = "FluidSystem_" + hashIt(typeName)
    generator = SimpleGenerator("FluidSystem", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName)
    return module.FluidSystem()
