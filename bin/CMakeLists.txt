# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

install(PROGRAMS
  testing/runtest.py
  RENAME dumux_runtest.py
  DESTINATION ${CMAKE_INSTALL_BINDIR})

install(PROGRAMS
  testing/fuzzycomparevtu.py
  testing/fuzzycomparedata.py
  DESTINATION ${CMAKE_INSTALL_BINDIR})
