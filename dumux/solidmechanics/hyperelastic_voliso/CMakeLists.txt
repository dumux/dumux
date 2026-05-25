# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_HEADERS *.hh *.inc)
install(FILES ${DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/solidmechanics/hyperelastic_voliso)
