from dune.typeregistry import generateTypeName
from dune.generator import Method
from dune.generator.generator import simpleGenerator

def OnePProblem(gridView, typetag):
    typeName, includes = generateTypeName("Dumux::OnePTestProblem", typetag)
    includes += ["dumux/porousmediumflow/py/onepproblem.hh", \
                 "dumux/porousmediumflow/py/problem.hh"]
    # temperature = Method("temperature", "&DuneType::temperature")
    return simpleGenerator([], "OnePProblem", "Dumux")(includes, typeName).OnePProblem()

def OnePTpfa(gridView):
    typeName, includes = generateTypeName("Dumux::TempTypeTag", gridView._typeName )
    includes += ["dumux/porousmediumflow/properties.hh", "dumux/porousmediumflow/py/oneptpfa.hh"]
    return simpleGenerator([], "OnePTpfa", "Dumux")(includes, typeName).OnePTpfa()
