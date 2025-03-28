# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

install(PROGRAMS dumux_runtest.py dumux_fuzzycompare_legacy.py
  DESTINATION ${CMAKE_INSTALL_BINDIR}
)
