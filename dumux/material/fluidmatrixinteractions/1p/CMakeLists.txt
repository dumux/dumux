# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_MATERIAL_FLUIDMATRIXINTERACTIONS_1P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MATERIAL_FLUIDMATRIXINTERACTIONS_1P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/fluidmatrixinteractions/1p)
