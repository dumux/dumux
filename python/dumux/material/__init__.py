from dune.generator.generator import SimpleGenerator

def Component(Type, Scalar="double"):
    moduleName = Type
    typeName = "Dumux::Components::" + Type + "<" + Scalar +  ">"
    includes = ["dumux/python/material/component.hh"]

    if Type == "Air":
        includes += ["dumux/material/components/air.hh"]
    elif Type == "Calcite":
        includes += ["dumux/material/components/calcite.hh"]
    elif Type == "H2O":
        includes += ["dumux/material/components/h2o.hh"]
    elif Type == "N2":
        includes += ["dumux/material/components/n2.hh"]

    generator = SimpleGenerator("Component", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName)

    return module.Component()
