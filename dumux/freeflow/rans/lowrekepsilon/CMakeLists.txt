# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_FREEFLOW_RANS_LOWREKEPSILON_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_RANS_LOWREKEPSILON_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans/lowrekepsilon)
