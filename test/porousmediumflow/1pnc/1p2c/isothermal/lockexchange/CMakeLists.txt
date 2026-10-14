# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_1p2c_lockexchange_boussinesq_box
              LABELS porousmediumflow 1p2c 1pnc
              SOURCES main.cc
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_lockexchange_boussinesq_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_boussinesq_box-00100.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_boussinesq_box params.input -Problem.Name test_1p2c_lockexchange_boussinesq_box")
