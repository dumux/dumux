# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "plot_results.py")

dumux_add_test_target(NAME test_2p1cni_heatpipe_box
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=TwoPOneCHeatPipeBox")

dumux_add_test_target(NAME test_2p1cni_heatpipe_tpfa
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=TwoPOneCHeatPipeTpfa")

dumux_add_test(NAME test_2p1cni_heatpipe_box
               TARGET test_2p1cni_heatpipe_box
               LABELS porousmediumflow 2p1c 2p1cni
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_heatpipe_box params.input -Problem.Name test_2p1cni_heatpipe_box"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_2p1cni_heatpipe_box-errors.csv
                        --check "L2(S) mean 0.8 1.2"
                        --check "L2(p_g) mean 0.8 1.2"
                        --check "L2(p_l) mean 0.8 1.2"
                        --check "L2(T) mean 0.8 1.2")

dumux_add_test(NAME test_2p1cni_heatpipe_tpfa
               TARGET test_2p1cni_heatpipe_tpfa
               LABELS porousmediumflow 2p1c 2p1cni
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_heatpipe_tpfa params.input -Problem.Name test_2p1cni_heatpipe_tpfa"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_2p1cni_heatpipe_tpfa-errors.csv
                        --check "L2(S) mean 0.8 1.2"
                        --check "L2(p_g) mean 0.8 1.2"
                        --check "L2(p_l) mean 0.8 1.2"
                        --check "L2(T) mean 0.8 1.2")

# central weighting of the mobilities makes both schemes second order
dumux_add_test(NAME test_2p1cni_heatpipe_box_central
               TARGET test_2p1cni_heatpipe_box
               LABELS porousmediumflow 2p1c 2p1cni
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_heatpipe_box params.input -Problem.Name test_2p1cni_heatpipe_box_central -Flux.UpwindWeight 0.5"
                        --levels 0 1 2 3 4 5
                        --level-args "-Grid.Refinement {}"
                        --errors test_2p1cni_heatpipe_box_central-errors.csv
                        --check "L2(S) mean 1.8 2.2"
                        --check "L2(p_g) mean 1.8 2.2"
                        --check "L2(p_l) mean 1.8 2.2"
                        --check "L2(T) mean 1.8 2.2")

dumux_add_test(NAME test_2p1cni_heatpipe_tpfa_central
               TARGET test_2p1cni_heatpipe_tpfa
               LABELS porousmediumflow 2p1c 2p1cni
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p1cni_heatpipe_tpfa params.input -Problem.Name test_2p1cni_heatpipe_tpfa_central -Flux.UpwindWeight 0.5"
                        --levels 0 1 2 3 4 5
                        --level-args "-Grid.Refinement {}"
                        --errors test_2p1cni_heatpipe_tpfa_central-errors.csv
                        --check "L2(S) mean 1.8 2.2"
                        --check "L2(p_g) mean 1.8 2.2"
                        --check "L2(p_l) mean 1.8 2.2"
                        --check "L2(T) mean 1.8 2.2")
