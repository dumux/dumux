# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

#install headers
install(FILES
biominsimplechemistry.hh
icpcomplexsalinitybrine.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/fluidsystems)
