# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_MULTIDOMAIN_FACET_CELLCENTERED_MPFA_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_FACET_CELLCENTERED_MPFA_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/facet/cellcentered/mpfa)
