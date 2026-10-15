# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

install(FILES
    model.hh
    localresidual.hh
    variables.hh
    couplingmanager.hh
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/poromechanics/poroelastic_large_deformations)
