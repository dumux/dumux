# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(expression)
add_subdirectory(format)
add_subdirectory(grid)
add_subdirectory(json)
add_subdirectory(vtk)
add_subdirectory(xml)

file(GLOB DUMUX_IO_HEADERS *.hh *.inc)
install(FILES ${DUMUX_IO_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/io)
