# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Classes and functions related to discretization methods"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


@cppWrapperCreator
def _createGridGeometry(gridView, discMethod="cctpfa"):
    """Construct a GridGeometry from a gridView"""

    includes = gridView._includes + ["dumux/python/discretization/gridgeometry.hh"]

    if discMethod == "cctpfa":
        includes += ["dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh"]
        typeName = "Dumux::CCTpfaFVGridGeometry<" + gridView._typeName + ">"
    elif discMethod == "box":
        includes += ["dumux/discretization/box/fvgridgeometry.hh"]
        typeName = "Dumux::BoxFVGridGeometry<double, " + gridView._typeName + ">"
    else:
        raise ValueError(f"Unknown discMethod {discMethod}")

    moduleName = "gridgeometry_" + hashIt(typeName)
    holderType = f"std::shared_ptr<{typeName}>"
    generator = SimpleGenerator("GridGeometry", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName, options=[holderType])
    return module.GridGeometry(gridView)


@cppWrapperClassAlias(creator=_createGridGeometry)
class GridGeometry:
    """Class alias used to instantiate a GridGeometry"""


@cppWrapperCreator
def _createGridVariables(*, problem, model):
    """Construct a GridGeometry from problem and the model"""

    includes = [
        "dumux/discretization/fvgridvariables.hh",
        "dumux/python/discretization/gridvariables.hh",
    ]
    ggeo = f"Dumux::GetPropType<{model.cppType}, Dumux::Properties::GridGeometry>"
    gvv = f"Dumux::GetPropType<{model.cppType}, Dumux::Properties::GridVolumeVariables>"
    gfc = f"Dumux::GetPropType<{model.cppType}, Dumux::Properties::GridFluxVariablesCache>"
    typeName = f"Dumux::FVGridVariables<{ggeo}, {gvv}, {gfc}>"

    moduleName = "gridvariables_" + hashIt(typeName)
    holderType = f"std::shared_ptr<{typeName}>"
    generator = SimpleGenerator("GridVariables", "Dumux::Python")
    module = generator.load(
        includes, typeName, moduleName, options=[holderType], preamble=model.cppHeader
    )
    module.GridVariables.model = property(lambda self: model)
    return module.GridVariables(problem, problem.gridGeometry())


@cppWrapperClassAlias(creator=_createGridVariables)
class GridVariables:
    """Class alias used to instantiate GridVariables"""
