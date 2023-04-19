# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_PORENETWORK_2P_STATIC_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PORENETWORK_2P_STATIC_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porenetwork/2p/static)
