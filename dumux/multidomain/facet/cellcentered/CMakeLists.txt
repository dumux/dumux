# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(mpfa)
add_subdirectory(tpfa)

file(GLOB DUMUX_MULTIDOMAIN_FACET_CELLCENTERED_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_FACET_CELLCENTERED_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/facet/cellcentered)
