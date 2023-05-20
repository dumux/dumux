# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES grids "params.input")

add_executable(test_porenetwork_solid_energy EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_porenetwork_solid_energy_stationary
              LABELS porenetwork solid energy
              TARGET test_porenetwork_solid_energy
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_solid_energy-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_solid_energy_stationary-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_porenetwork_solid_energy params.input
                                    -Problem.Name test_pnm_solid_energy_stationary -Problem.IsStationary true")

dumux_add_test(NAME test_porenetwork_solid_energy_transient
              LABELS porenetwork solid energy
              TARGET test_porenetwork_solid_energy
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_solid_energy-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_solid_energy_transient-00033.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_porenetwork_solid_energy params.input
                                    -Problem.Name test_pnm_solid_energy_transient -Problem.IsStationary false")
