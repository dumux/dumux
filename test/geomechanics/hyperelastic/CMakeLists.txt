# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")
dumux_add_test(
    NAME test_hyperelastic_box
    LABELS geomechanics elastic
    SOURCES main.cc
    COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_box-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_box-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_box params.input"
)
