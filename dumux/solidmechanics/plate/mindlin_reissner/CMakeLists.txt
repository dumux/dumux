# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_SOLIDMECHANICS_PLATE_MINDLIN_REISSNER_HEADERS *.hh)
install(FILES ${DUMUX_SOLIDMECHANICS_PLATE_MINDLIN_REISSNER_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/solidmechanics/plate/mindlin_reissner)
