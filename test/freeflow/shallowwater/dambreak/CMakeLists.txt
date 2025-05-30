# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

add_executable(test_shallowwater_dambreak EXCLUDE_FROM_ALL main.cc)

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
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_dambreak-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0002-dambreak_parallel-00001.pvtu
                             --zeroThreshold {"velocityY":1e-14}
                             --ignore "process rank"
                             --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_dambreak params.input
                           -Problem.Name dambreak_parallel")
