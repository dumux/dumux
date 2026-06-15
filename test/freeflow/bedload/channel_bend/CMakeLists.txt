# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES grid "params.input")

dumux_add_test(NAME test_bedload_channel_bend
               SOURCES main.cc
               LABELS bedload
               TIMEOUT 40
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS   --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_bedload_channel_bend.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/channel_bend-00006.vtu
                          --zeroThreshold {"velocityY":1e-8,"layerThickness":1e-12}
                          --command "${CMAKE_CURRENT_BINARY_DIR}/test_bedload_channel_bend params.input")
