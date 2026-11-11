# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

#add links to input files
dune_symlink_to_source_files(FILES "params_box.input" "params_tpfa.input")

dumux_add_test_target(NAME test_2p1cni_steaminjection_tpfa
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=TwoPOneCNICCTpfa")

dumux_add_test_target(NAME test_2p1cni_steaminjection_box
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=TwoPOneCNIBox")

dumux_add_test(NAME test_2p1cni_steaminjection_waterwet_box
              LABELS porousmediumflow 2p1c  2p1cni
              TARGET test_2p1cni_steaminjection_box
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p1cni_steaminjection_waterwet_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_waterwet_box-00007.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_box params_box.input -Problem.Name test_2p1cni_steaminjection_waterwet_box")

dumux_add_test(NAME test_2p1cni_steaminjection_waterwet_tpfa
              LABELS porousmediumflow 2p1c  2p1cni
              TARGET test_2p1cni_steaminjection_tpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p1cni_steaminjection_waterwet_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_waterwet_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_tpfa params_tpfa.input -Problem.Name test_2p1cni_steaminjection_waterwet_tpfa")

dumux_add_test(NAME test_2p1cni_steaminjection_gaswet_box
              LABELS porousmediumflow 2p1c  2p1cni
              TARGET test_2p1cni_steaminjection_box
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p1cni_steaminjection_gaswet_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_gaswet_box-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_box params_box.input
                                                                               -Problem.Name test_2p1cni_steaminjection_gaswet_box
                                                                               -SpatialParams.GasWetting true")

dumux_add_test(NAME test_2p1cni_steaminjection_gaswet_tpfa
              LABELS porousmediumflow 2p1c  2p1cni
              TARGET test_2p1cni_steaminjection_tpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p1cni_steaminjection_gaswet_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_gaswet_tpfa-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_steaminjection_tpfa params_tpfa.input
                                                                                -Problem.Name test_2p1cni_steaminjection_gaswet_tpfa
                                                                                -SpatialParams.GasWetting true")
