# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(boundary)
add_subdirectory(dualnetwork)
add_subdirectory(embedded)
add_subdirectory(facet)
add_subdirectory(freeflow)
add_subdirectory(io)

file(GLOB DUMUX_MULTIDOMAIN_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain)
