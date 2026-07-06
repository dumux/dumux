# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_PHASEFIELD_CAHNHILLIARD_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PHASEFIELD_CAHNHILLIARD_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/phasefield/cahnhilliard)
