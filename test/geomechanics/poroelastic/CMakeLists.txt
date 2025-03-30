# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# using box and numeric differentiation
dumux_add_test(NAME test_poroelastic_box
              LABELS geomechanics poroelastic
              SOURCES main.cc
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_poroelastic_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_box params.input -Problem.Name test_poroelastic_box")

dumux_add_test(NAME test_poroelastic_box_gridformat_io
               LABELS geomechanics poroelastic
               CMAKE_GUARD DUMUX_HAVE_GRIDFORMAT
               COMPILE_DEFINITIONS USE_GRIDFORMAT=1
               SOURCES main.cc
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS  --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_poroelastic_box-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_box_gridformat_io-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_box_gridformat_io params.input -Problem.Name test_poroelastic_box_gridformat_io")
