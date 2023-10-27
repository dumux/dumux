# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(1p)
add_subdirectory(1pnc)
add_subdirectory(2p)
add_subdirectory(common)
add_subdirectory(solidenergy)
add_subdirectory(util)

file(GLOB DUMUX_PORENETWORK_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PORENETWORK_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porenetwork)
