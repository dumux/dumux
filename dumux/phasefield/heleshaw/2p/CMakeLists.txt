# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_PHASEFIELD_HELESHAW_2P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PHASEFIELD_HELESHAW_2P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/phasefield/heleshaw/2p)

add_subdirectory(cvfe)
