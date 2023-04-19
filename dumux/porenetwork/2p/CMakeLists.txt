# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(static)

file(GLOB DUMUX_PORENETWORK_2P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PORENETWORK_2P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porenetwork/2p)
