# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_MATERIAL_COMPONENTS_IAPWS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MATERIAL_COMPONENTS_IAPWS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/components/iapws)
