# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(limitednikuradse)
dune_symlink_to_source_files(FILES params.input convergencetest.py plot_results.py)

dumux_add_test(NAME test_shallowwater_roughchannel
               SOURCES main.cc
               LABELS shallowwater
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_roughchannel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/roughchannel-00001.vtu
                             --zeroThreshold {"velocityY":1e-14}
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_roughchannel params.input
                             -Problem.Name roughchannel")

# convergence towards the normal flow
dumux_add_test(NAME test_shallowwater_roughchannel_convergence
               TARGET test_shallowwater_roughchannel
               LABELS shallowwater
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS test_shallowwater_roughchannel)
