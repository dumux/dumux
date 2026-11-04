# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(momentum)

dune_symlink_to_source_files(FILES params.input)
dumux_add_test(NAME test_ff_stokes_channel_pipe
               LABELS freeflow navierstokes
               SOURCES main.cc
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pipe params.input -Problem.Name test_ff_stokes_channel_pipe -Problem.PrintConvergenceTestFile true"
                        --levels 0 1 2
                        --level-args "-Grid.Refinement {}"
                        --errors test_ff_stokes_channel_pipe-errors.csv
                        --check "L2Abs(p) mean 1.8"
                        --check "L2Abs(vx) mean 1.9"
                        --check "L2Abs(vy) mean 1.9")
