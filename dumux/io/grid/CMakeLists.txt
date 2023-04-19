# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(porenetwork)

file(GLOB DUMUX_IO_GRID_HEADERS *.hh *.inc)
install(FILES ${DUMUX_IO_GRID_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/io/grid)
