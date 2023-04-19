# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(thermal)

file(GLOB DUMUX_POROUSMEDIUMFLOW_NONEQUILIBRIUM_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_NONEQUILIBRIUM_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/nonequilibrium)
