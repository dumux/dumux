# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(binarycoefficients)
add_subdirectory(chemistry)
add_subdirectory(components)
add_subdirectory(constraintsolvers)
add_subdirectory(eos)
add_subdirectory(fluidmatrixinteractions)
add_subdirectory(fluidstates)
add_subdirectory(fluidsystems)
add_subdirectory(solidstates)
add_subdirectory(solidsystems)

file(GLOB DUMUX_MATERIAL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MATERIAL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material)
