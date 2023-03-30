# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Components are building blocks of fluid and solid systems"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


_components = {
    "Air": "dumux/material/components/air.hh",
    "H2O": "dumux/material/components/h2o.hh",
    "SimpleH2O": "dumux/material/components/simpleh2o.hh",
    "N2": "dumux/material/components/n2.hh",
    "Calcite": "dumux/material/components/calcite.hh",
    "Constant": "dumux/material/components/constant.hh",
}


def listComponents():
    """List all available component names"""
    print("The following components are available:")
    print(_components.keys())


@cppWrapperCreator
def _createComponent(name, *, scalar="double", componentId=0):
    """Create a new component of the given name"""
    if name == "Constant":
        typeName = f"Dumux::Components::{name} <{componentId}, {scalar}>"
    else:
        typeName = f"Dumux::Components::{name} <{scalar}>"
    moduleName = f"{name.lower()}_{hashIt(typeName)}"
    includes = ["dumux/python/material/components/component.hh"]
    includes += [_components[name]]
    generator = SimpleGenerator("Component", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName)
    return module.Component()


@cppWrapperClassAlias(creator=_createComponent)
class Component:
    """Class alias used to instantiate a Component"""
