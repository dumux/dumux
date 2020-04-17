from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

# construct a GridGeometry from a gridView
# the grid geometry is JIT compiled
def GridGeometry(gridView, discMethod="cctpfa"):
    includes = gridView._includes + ["dumux/python/discretization/gridgeometry.hh"]

    if discMethod == "cctpfa":
        includes += ["dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh"]
        typeName = "Dumux::CCTpfaFVGridGeometry<" + gridView._typeName + ">"
    elif discMethod == "box":
        includes += ["dumux/discretization/box/fvgridgeometry.hh"]
        typeName = "Dumux::BoxFVGridGeometry<double, " + gridView._typeName + ">"
    else:
        raise ValueError("Unknown discMethod {}".format(discMethod))

    moduleName = "gridgeometry_" + hashIt(typeName)
    holderType = "std::shared_ptr<{}>".format(typeName)
    generator = SimpleGenerator("GridGeometry", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName, options=[holderType])
    return module.GridGeometry(gridView)
