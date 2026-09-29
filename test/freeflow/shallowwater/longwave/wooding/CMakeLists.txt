# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_vegetated.input)

add_executable(test_shallowwater_longwave_wooding EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_shallowwater_longwave_wooding PUBLIC TYPETAG=WoodingLongWave)

# regression tests for the three long-wave approximations and the local-inertial shallow water equations
foreach(approximation diffusive kinematic inertiacorrected)
    dumux_add_test(NAME test_shallowwater_longwave_wooding_${approximation}
                   TARGET test_shallowwater_longwave_wooding
                   LABELS shallowwater
                   CMAKE_GUARD HAVE_UMFPACK
                   COMMAND ${DUMUX_RUNTEST}
                   CMD_ARGS --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_longwave_wooding_${approximation}-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/wooding_${approximation}-00001.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_longwave_wooding params.input
                                       -LongWave.Approximation ${approximation} -Problem.Name wooding_${approximation}")
endforeach()

# the depth is probed on the second row of cells, since the first resolves the drawdown of the free overfall
dumux_add_test(NAME test_shallowwater_localinertial_wooding
               SOURCES main.cc
               LABELS shallowwater
               COMPILE_DEFINITIONS TYPETAG=WoodingLocalInertial
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_localinertial_wooding-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/wooding_localinertial-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_localinertial_wooding params.input
                                   -Problem.Name wooding_localinertial
                                   -Problem.DepthProbeMinY 25 -Problem.DepthProbeMaxY 50
                                   -Assembly.NumericDifference.PriVarMagnitude \"1 1 1\"")

# the vegetated catchment reaches equilibrium, where the outflow equals the rainfall on the catchment
foreach(approximation diffusive inertiacorrected)
    dumux_add_test(NAME test_shallowwater_longwave_wooding_vegetated_${approximation}
                   TARGET test_shallowwater_longwave_wooding
                   LABELS shallowwater
                   CMAKE_GUARD HAVE_UMFPACK
                   CMD_ARGS params_vegetated.input
                            -LongWave.Approximation ${approximation}
                            -Problem.Name wooding_vegetated_${approximation}
                            -Problem.EquilibriumTolerance 1e-6)
endforeach()

# The kinematic wave carries no flux across the channel, whose bed is level across. Water
# arriving at the outlet in the channel's corners leaves only through their share of the
# outlet, so the corners fill slowly and the outflow is still below equilibrium when the rain stops.
dumux_add_test(NAME test_shallowwater_longwave_wooding_vegetated_kinematic
               TARGET test_shallowwater_longwave_wooding
               LABELS shallowwater
               CMAKE_GUARD HAVE_UMFPACK
               CMD_ARGS params_vegetated.input
                        -LongWave.Approximation kinematic
                        -Problem.Name wooding_vegetated_kinematic)

dumux_add_test(NAME test_shallowwater_localinertial_wooding_vegetated
               TARGET test_shallowwater_localinertial_wooding
               LABELS shallowwater
               CMAKE_GUARD HAVE_UMFPACK
               CMD_ARGS params_vegetated.input
                        -Problem.Name wooding_vegetated_localinertial
                        -Problem.DepthProbeMinY 25 -Problem.DepthProbeMaxY 50
                        -Assembly.NumericDifference.PriVarMagnitude "1 1 1"
                        -Problem.EquilibriumTolerance 1e-6)
