# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(box)
add_subdirectory(ccmpfa)
add_subdirectory(cctpfa)
add_subdirectory(cvfe)
add_subdirectory(porenetwork)
add_subdirectory(shallowwater)
add_subdirectory(staggered)

file(GLOB DUMUX_FLUX_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FLUX_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/flux)
