# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB DUMUX_PYTHON_IO_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PYTHON_IO_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/python/io)
