# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# the exact solution is known, so the tests check themselves instead of comparing with a reference

dumux_add_test(NAME test_shallowwater_lakeatrest
               SOURCES main.cc
               LABELS shallowwater
               CMD_ARGS params.input
                        -Problem.ErrorTolerance 1e-10)

# on a partially dry bed the accuracy is limited by the thin film in the dry cells,
# but still orders of magnitude below the bed variation a non-well-balanced scheme would produce
dumux_add_test(NAME test_shallowwater_lakeatrest_emerged
               TARGET test_shallowwater_lakeatrest
               LABELS shallowwater
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_lakeatrest
               CMD_ARGS params.input
                        -Problem.Name lakeatrest-emerged
                        -Problem.WaterLevel 0.5
                        -Problem.MinWaterDepth 1e-5
                        -Problem.ErrorTolerance 1e-5)

dumux_add_test(NAME test_shallowwater_lakeatrest_parallel
               TARGET test_shallowwater_lakeatrest
               LABELS shallowwater
               CMAKE_GUARD MPI_FOUND
               COMMAND ${MPIEXEC}
               CMD_ARGS -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_lakeatrest
                        params.input
                        -Problem.Name lakeatrest-parallel
                        -Problem.ErrorTolerance 1e-10)
