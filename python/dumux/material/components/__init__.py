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
}


def listComponents():
    """List all available component names"""
    print("The following components are availabe:")
    print(_components.keys())


@cppWrapperCreator
def _createComponent(name, *, scalar="double"):
    """Create a new component of the given name"""

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
