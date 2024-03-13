# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_MULTIDOMAIN_PORENETWORK_CONSTRAINT_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_PORENETWORK_CONSTRAINT_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/porenetwork/constraint)
