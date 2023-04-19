# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(box)
add_subdirectory(cellcentered)
add_subdirectory(cvfe)
add_subdirectory(facecentered)
add_subdirectory(fem)
add_subdirectory(porenetwork)
add_subdirectory(pq1bubble)
add_subdirectory(projection)
add_subdirectory(staggered)

file(GLOB DUMUX_DISCRETIZATION_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DISCRETIZATION_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization)
