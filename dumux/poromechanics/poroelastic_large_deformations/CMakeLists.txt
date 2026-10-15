# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/poromechanics/poroelastic_large_deformations)
