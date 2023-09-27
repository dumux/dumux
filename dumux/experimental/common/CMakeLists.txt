# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(typetraits)
file(GLOB DUMUX_EXPERIMENTAL_COMMON_HEADERS *.hh *.inc)
install(FILES ${DUMUX_EXPERIMENTAL_COMMON_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/experimental/common)
