# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_rans_zeroeq_channel
              LABELS freeflow navierstokes rans
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              TIMEOUT 600
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_rans_zeroeq_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_rans_zeroeq_channel-00046.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_rans_zeroeq_channel params.input")

add_executable(test_rans_zeroeq_channelni EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_rans_zeroeq_channelni PUBLIC "NONISOTHERMAL=1")

dumux_add_test(NAME test_rans_zeroeq_channelni
              TARGET test_rans_zeroeq_channelni
              LABELS freeflow navierstokes rans
              CMAKE_GUARD HAVE_UMFPACK
              TIMEOUT 600
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_rans_zeroeq_channelni-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_rans_zeroeq_channelni-00053.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_rans_zeroeq_channelni params_ni.input")
