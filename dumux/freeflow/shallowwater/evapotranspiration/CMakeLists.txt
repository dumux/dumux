# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_FREEFLOW_SHALLOWWATER_EVAPOTRANSPIRATION_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_SHALLOWWATER_EVAPOTRANSPIRATION_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/shallowwater/evapotranspiration)
