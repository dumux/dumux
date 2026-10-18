# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params_2d.input" "params_3d.input")

dumux_add_test(NAME test_2pve_immiscible_tpfa_2d
               SOURCES main.cc
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
               LABELS porousmediumflow 2pve
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --zeroThreshold "{\"S_gas\":2e-3,\"mob_gas\":90,\"velocity_gas (m/s)\":1e-7}"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2pve_immiscible_tpfa_2d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/2pve_lenses_2d_40x0x18-00010.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pve_immiscible_tpfa_2d params_2d.input -Problem.Name 2pve_lenses_2d_ -Grid.LowerLeft '0 0' -Grid.UpperRight '250 55' -Grid.Cells '40 18' -TimeLoop.TEnd 4e5 -TimeLoop.MaxTimeStepSize 10000 -BoundaryConditions.InjectionRate -5.835e-4")

# Homogeneous formation, injection stops at 8e4 s: 40 columns, fine cells of 1 m height.
# Prescribed startup ramp and fixed steps keep trapping history independent of Newton iteration counts.
# Require at least 1% of the injected gas to remain residually trapped.
dumux_add_test(NAME test_2pve_immiscible_tpfa_2d_trapping
               TARGET test_2pve_immiscible_tpfa_2d
               LABELS porousmediumflow 2pve
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --zeroThreshold "{\"S_gas\":3e-3,\"mob_gas\":90,\"velocity_gas (m/s)\":1e-7,\"velocity_liq (m/s)\":1e-7}"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2pve_immiscible_tpfa_2d_trapping_coarse-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/2pve_trapping_2d_40x0x55-00011.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_2pve_immiscible_tpfa_2d_trapping_fine-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/fine_2pve_trapping_2d_40x0x55-00011.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pve_immiscible_tpfa_2d params_2d.input -Problem.Name 2pve_trapping_2d_ -Grid.LowerLeft '0 0' -Grid.UpperRight '250 55' -Grid.Cells '40 55' -TimeLoop.TEnd 2e6 -TimeLoop.MaxTimeStepSize 390 -TimeLoop.UsePrescribedTimeSteps true -BoundaryConditions.InjectionRate -5.835e-4 -BoundaryConditions.InjectionEndTime 8e4 -SpatialParams.Snr 0.2 -SpatialParams.EnableLens false -MassBalance.MinTrappedGasFraction 0.01")

dumux_add_test(NAME test_2pve_immiscible_tpfa_3d
               SOURCES main.cc
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>
               LABELS porousmediumflow 2pve
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --zeroThreshold "{\"S_gas\":6e-4,\"mob_gas\":10,\"velocity_gas (m/s)\":1e-8}"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2pve_immiscible_tpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/2pve_lenses_3d_15x8x3-00010.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pve_immiscible_tpfa_3d params_3d.input -Problem.Name 2pve_lenses_3d_ -Grid.LowerLeft '0 0 0' -Grid.UpperRight '250 100 25' -Grid.Cells '15 8 3' -TimeLoop.TEnd 2e6 -TimeLoop.MaxTimeStepSize 10000 -BoundaryConditions.InjectionRate -1.45875e-3")

# the well lies at the center of the y-extent also if the domain does not start at y = 0
dumux_add_test(NAME test_2pve_immiscible_tpfa_3d_shifted
               TARGET test_2pve_immiscible_tpfa_3d
               LABELS porousmediumflow 2pve
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_2pve_immiscible_tpfa_3d
               CMD_ARGS params_3d.input
                        -Problem.Name 2pve_shifted_3d_
                        -Grid.LowerLeft "0 100 0"
                        -Grid.UpperRight "250 200 25"
                        -TimeLoop.TEnd 0.125)
