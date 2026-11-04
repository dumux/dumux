# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)
dumux_add_test(NAME test_ff_stokes_channel_pipe_momentum_diamond
               LABELS freeflow navierstokes
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=PipeFlowDiamond
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pipe_momentum_diamond params.input -Problem.Name test_ff_stokes_channel_pipe_momentum_diamond"
                        --levels 0 1 2
                        --level-args "-Grid.Refinement {}"
                        --errors test_ff_stokes_channel_pipe_momentum_diamond-errors.csv
                        --check "L2Abs(vx) mean 1.9"
                        --check "L2Abs(vy) mean 1.9")

dumux_add_test(NAME test_ff_stokes_channel_pipe_momentum_pq1bubble
               LABELS freeflow navierstokes
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=PipeFlowPQ1Bubble
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pipe_momentum_pq1bubble params.input -Problem.Name test_ff_stokes_channel_pipe_momentum_pq1bubble"
                        --levels 0 1 2
                        --level-args "-Grid.Refinement {}"
                        --errors test_ff_stokes_channel_pipe_momentum_pq1bubble-errors.csv
                        --check "L2Abs(vx) mean 1.9")
