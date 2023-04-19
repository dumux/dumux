# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES params.input)

dumux_add_test(NAME test_2p3c_surfactant
               LABELS porousmediumflow 2pnc
               SOURCES main.cc
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p3c_surfactant-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p3c_surfactant-00144.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p3c_surfactant params.input -Problem.Name test_2p3c_surfactant")
