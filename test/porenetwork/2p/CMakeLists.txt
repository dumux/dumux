# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(static)

add_input_file_links()

dumux_add_test(NAME test_pnm_2p
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMMAND ${DUMUX_RUNTEST}
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00093.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p")

dumux_add_test(NAME test_pnm_2p_without_corner_flow
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMMAND ${DUMUX_RUNTEST}
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p_without_corner_flow-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_without_corner_flow-00114.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_without_corner_flow params.input
                                        -Problem.Name test_pnm_2p_without_corner_flow -SpatialParams.ContactAngle 0.8")

dumux_add_test(NAME test_pnm_2p_accuracy_check
               TARGET test_pnm_2p
               LABELS porenetwork
               COMMAND ${DUMUX_RUNTEST}
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p_accuracy_check-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_accuracy_check-00171.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p params.input
                                        -Problem.Name test_pnm_2p_accuracy_check
                                        -TimeLoop.MaxTimeStepSize 1e-4 -TimeLoop.TEnd 0.005
                                        -InvasionState.AccuracyCriterion 0.99")

dumux_add_test(NAME test_pnm_2p_reg
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 USETHETAREGULARIZATION=1
               COMMAND ${DUMUX_RUNTEST}
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p_reg-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_reg-00259.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_reg params.input
                                        -Problem.Name test_pnm_2p_reg
                                        -TimeLoop.MaxTimeStepSize 1e-4 -TimeLoop.TEnd 0.005
                                        -InvasionState.AccuracyCriterion 0.99
                                        -Problem.RegularizationDelta 1e-4
                                        -InvasionState.InvasionThetaThreshold 1e-14")

dumux_add_test(NAME test_pnm_2pni
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0
               COMMAND ${DUMUX_RUNTEST}
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni-00133.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni params_ni.input -Problem.Name test_pnm_2pni")
