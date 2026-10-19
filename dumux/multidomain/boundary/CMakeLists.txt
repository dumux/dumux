# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(darcydarcy)
add_subdirectory(freeflowporenetwork)
add_subdirectory(freeflowporousmedium)

file(GLOB DUMUX_MULTIDOMAIN_BOUNDARY_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_BOUNDARY_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/boundary)
