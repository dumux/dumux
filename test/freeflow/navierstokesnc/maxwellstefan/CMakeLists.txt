
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_ff_stokes2c_maxwellstefan
              LABELS freeflow navierstokesnc navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_maxwellstefan-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_maxwellstefan-00005.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_maxwellstefan params.input
                             -Problem.Name test_ff_stokes2c_maxwellstefan
                             -Problem.PlotOutput false")
