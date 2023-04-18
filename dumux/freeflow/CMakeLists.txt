# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(compositional)
add_subdirectory(navierstokes)
add_subdirectory(nonisothermal)
add_subdirectory(rans)
add_subdirectory(shallowwater)

file(GLOB DUMUX_FREEFLOW_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow)
