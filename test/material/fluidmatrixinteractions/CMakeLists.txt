# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(2p)
add_subdirectory(porenetwork)

dumux_add_test(SOURCES test_effectivediffusivity.cc
              LABELS unit material
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzyData --delimiter " "
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_effectivediffusivityconstanttau-reference.dat
                               ${CMAKE_CURRENT_BINARY_DIR}/constant_d_eff.dat
                               ${CMAKE_SOURCE_DIR}/test/references/test_effectivediffusivitymillingtonquirk-reference.dat
                               ${CMAKE_CURRENT_BINARY_DIR}/millingtonquirk_d_eff.dat
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_effectivediffusivity -SpatialParams.Tortuosity 0.3")
