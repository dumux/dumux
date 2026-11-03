# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_1p_internaldirichlet_tpfa
              SOURCES main.cc
              LABELS porousmediumflow 1p
              COMPILE_DEFINITIONS TYPETAG=OnePInternalDirichletTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_internaldirichlet_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_internaldirichlet_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_internaldirichlet_tpfa params.input -Problem.Name test_1p_internaldirichlet_tpfa")

dumux_add_test(NAME test_1p_internaldirichlet_box
              SOURCES main.cc
              LABELS porousmediumflow 1p
              COMPILE_DEFINITIONS TYPETAG=OnePInternalDirichletBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_internaldirichlet_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_internaldirichlet_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_internaldirichlet_box params.input -Problem.Name test_1p_internaldirichlet_box")
