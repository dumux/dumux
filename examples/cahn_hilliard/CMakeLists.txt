# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# A simple test to automatically check the example compiles and produces the same result
dumux_add_test(NAME example_cahn_hilliard
               LABELS example
               SOURCES main.cc
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_cahn_hilliard-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cahn_hilliard-00102.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/example_cahn_hilliard params.input
                        -Grid.Cells \" 10 10\" -Problem.Mobility 0.01 -Problem.SurfaceTension 0.5")

dumux_add_test(NAME example_cahn_hilliard_parallel
               LABELS example
               TARGET example_cahn_hilliard
               CMAKE_GUARD MPI_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_cahn_hilliard_parallel_p0-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-cahn_hilliard_parallel-00102.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/example_cahn_hilliard_parallel_p1-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-cahn_hilliard_parallel-00102.vtu
                        --command "${MPIEXEC} -np 2
                        ${CMAKE_CURRENT_BINARY_DIR}/example_cahn_hilliard params.input
                        -Problem.Name cahn_hilliard_parallel
                        -Grid.Cells \" 10 10\" -Problem.Mobility 0.01 -Problem.SurfaceTension 0.5")
