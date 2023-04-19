# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_DISCRETIZATION_CELLCENTERED_MPFA_OMETHOD_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DISCRETIZATION_CELLCENTERED_MPFA_OMETHOD_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/cellcentered/mpfa/omethod)
