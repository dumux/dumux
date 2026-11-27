# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(common)
add_subdirectory(komega)
add_subdirectory(oneeq)
add_subdirectory(zeroeq)

file(GLOB DUMUX_FREEFLOW_RANS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_RANS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans)
