# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_shallowwater_roughchannel_limitednikuradse
               SOURCES main.cc
               LABELS shallowwater
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_roughchannel_limitednikuradse-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/roughchannel-00001.vtu
                             --zeroThreshold {"velocityY":1e-14}
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_roughchannel_limitednikuradse params.input
                                        -Problem.Name roughchannel")
