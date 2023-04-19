# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# To add an example it would be enough to have
# dumux_add_test(NAME example_diffusion SOURCES main.cc)
# But we also make sure with automated testing that the example keep compiling
# and produce the same result, so here is a regression test:
dumux_add_test(NAME example_diffusion
               LABELS example
               SOURCES main.cc
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_diffusion-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/diffusion-00050.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/example_diffusion params.input
                                   -Grid.Cells \"6 6\" -Problem.DiffusionCoefficient 0.001")

# We also add a parallel test to make sure the parallel version also works
dumux_add_test(NAME example_diffusion_parallel
               LABELS example
               CMAKE_GUARD MPI_FOUND
               TARGET example_diffusion
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_diffusion_parallel_p0-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-diffusion_parallel-00050.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/example_diffusion_parallel_p1-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-diffusion_parallel-00050.vtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/example_diffusion params.input
                                   -Grid.Cells \"6 6\" -Problem.DiffusionCoefficient 0.001 -Problem.Name diffusion_parallel")
