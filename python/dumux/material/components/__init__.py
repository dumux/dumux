from dune.generator.generator import SimpleGenerator


components = {
    "Air": "dumux/material/components/air.hh",
    "H2O": "dumux/material/components/h2o.hh",
    "SimpleH2O": "dumux/material/components/simpleh2o.hh",
    "N2": "dumux/material/components/n2.hh",
    "Calcite": "dumux/material/components/calcite.hh",
}


def Component(*, name, scalar="double"):
    moduleName = name
    typeName = f"Dumux::Components::{name} <{scalar}>"
    includes = ["dumux/python/material/components/component.hh"]
    includes += [components[name]]
    generator = SimpleGenerator("Component", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName)
    return module.Component()


def listComponents():
    print("The following components are availabe:")
    print(components.keys())
