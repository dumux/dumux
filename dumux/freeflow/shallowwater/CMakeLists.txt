# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(evapotranspiration)
add_subdirectory(longwave)
add_subdirectory(surfacerunoff)

file(GLOB DUMUX_FREEFLOW_SHALLOWWATER_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_SHALLOWWATER_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/shallowwater)
