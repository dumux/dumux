# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py" "plot_results.py")

add_executable(test_2p1cni_heatpipe_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p1cni_heatpipe_box PUBLIC "TYPETAG=TwoPOneCHeatPipeBox")

add_executable(test_2p1cni_heatpipe_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p1cni_heatpipe_tpfa PUBLIC "TYPETAG=TwoPOneCHeatPipeTpfa")

dumux_add_test(NAME test_2p1cni_heatpipe_box
               TARGET test_2p1cni_heatpipe_box
               LABELS porousmediumflow 2p1c 2p1cni
               COMMAND ./convergencetest.py
               CMD_ARGS test_2p1cni_heatpipe_box params.input)

dumux_add_test(NAME test_2p1cni_heatpipe_tpfa
               TARGET test_2p1cni_heatpipe_tpfa
               LABELS porousmediumflow 2p1c 2p1cni
               COMMAND ./convergencetest.py
               CMD_ARGS test_2p1cni_heatpipe_tpfa params.input)

# central weighting of the mobilities makes both schemes second order
dumux_add_test(NAME test_2p1cni_heatpipe_box_central
               TARGET test_2p1cni_heatpipe_box
               LABELS porousmediumflow 2p1c 2p1cni
               COMMAND ./convergencetest.py
               CMD_ARGS test_2p1cni_heatpipe_box params.input
                        --expected-rate 2 --refinements 6 --name test_2p1cni_heatpipe_box_central
                        -Flux.UpwindWeight 0.5)

dumux_add_test(NAME test_2p1cni_heatpipe_tpfa_central
               TARGET test_2p1cni_heatpipe_tpfa
               LABELS porousmediumflow 2p1c 2p1cni
               COMMAND ./convergencetest.py
               CMD_ARGS test_2p1cni_heatpipe_tpfa params.input
                        --expected-rate 2 --refinements 6 --name test_2p1cni_heatpipe_tpfa_central
                        -Flux.UpwindWeight 0.5)
