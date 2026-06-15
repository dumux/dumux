# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_bedload_bump
               SOURCES main.cc
               LABELS bedload
               TIMEOUT 40
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS   --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_bedload_bump.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/bump-00005.vtu
                          --zeroThreshold {"bedloadDischargeYGrainClass1":1e-8,"velocityY":1e-8,"layerThickness":1e-12}
                          --command "${CMAKE_CURRENT_BINARY_DIR}/test_bedload_bump params.input")
