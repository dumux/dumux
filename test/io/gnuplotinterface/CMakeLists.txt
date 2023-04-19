# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
dumux_add_test(SOURCES test_gnuplotinterface.cc
              LABELS unit io)
