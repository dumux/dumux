# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input plot_results.py)

dumux_add_test_target(NAME test_shallowwater_dambreak
                      SOURCES main.cc)

dumux_add_test(NAME test_shallowwater_dambreak_sequential
               TARGET test_shallowwater_dambreak
               LABELS shallowwater
               TIMEOUT 1500
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_dambreak-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/dambreak-00001.vtu
                             --zeroThreshold {"velocityY":1e-14}
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_dambreak params.input
                             -Problem.Name dambreak")

dumux_add_test(NAME test_shallowwater_dambreak_parallel
               TARGET test_shallowwater_dambreak
               LABELS shallowwater parallel
               TIMEOUT 1500
               CMAKE_GUARD MPI_FOUND
               MPI_RANKS 2
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_dambreak-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0002-dambreak_parallel-00001.pvtu
                             --zeroThreshold {"velocityY":1e-14}
                             --ignore "process rank"
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_dambreak params.input
                           -Problem.Name dambreak_parallel")

# convergence towards the exact solution of the Riemann problem, with the time step size proportional to the cell size
dumux_add_test(NAME test_shallowwater_dambreak_convergence
               TARGET test_shallowwater_dambreak
               LABELS shallowwater
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_dambreak params.input -Problem.Name test_shallowwater_dambreak_convergence -Grid.Cells1 1"
                        --levels "100 0.005" "200 0.0025" "400 0.00125"
                        --level-args "-Grid.Cells0 {0} -TimeLoop.MaxTimeStepSize {1} -TimeLoop.DtInitial {1}"
                        --errors test_shallowwater_dambreak_convergence-errors.csv
                        --check "L1Rel(depth) mean 0.7"
                        --check "L1Rel(u) mean 0.7")
