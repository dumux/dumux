# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""DuMux input-output library"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias


@cppWrapperCreator
def _createVtkOutputModule(*, gridVariables, solutionVector, name):
    """Construct a VtkOutputModule"""

    includes = gridVariables._includes + solutionVector._includes
    includes += ["dumux/python/io/vtkoutputmodule.hh", "dumux/io/vtkoutputmodule.hh"]
    typeName = f"Dumux::VtkOutputModule<{gridVariables._typeName}, {solutionVector._typeName}>"
    moduleName = "vtkoutputmodule_" + hashIt(typeName)
    generator = SimpleGenerator("VtkOutputModule", "Dumux::Python")
    module = generator.load(includes, typeName, moduleName, preamble=gridVariables.model.cppHeader)
    return module.VtkOutputModule(gridVariables, solutionVector, name)


@cppWrapperClassAlias(creator=_createVtkOutputModule)
class VtkOutputModule:
    """Class alias used to instantiate a VtkOutputModule"""
