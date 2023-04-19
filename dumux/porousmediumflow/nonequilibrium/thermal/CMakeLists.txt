# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_POROUSMEDIUMFLOW_NONEQUILIBRIUM_THERMAL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_NONEQUILIBRIUM_THERMAL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/nonequilibrium/thermal)
