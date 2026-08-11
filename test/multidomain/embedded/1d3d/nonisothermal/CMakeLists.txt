# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_wellbore_heat_transport
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/dumux_runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport -TimeLoop.TEnd 7200")


add_input_file_links()
dune_symlink_to_source_files(FILES "grids" "analytical_solution.py" "analytical_solution_detailed.py" "convergence_test.py")
