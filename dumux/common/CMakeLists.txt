# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(properties)
add_subdirectory(typetraits)

dune_library_add_sources(dumux SOURCES parameters.cc)

file(GLOB DUMUX_COMMON_HEADERS *.hh *.inc)
install(FILES ${DUMUX_COMMON_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/common)
