# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


dumux_add_test(NAME test_1p_internaldirichlet_tpfa
              SOURCES ${CMAKE_SOURCE_DIR}/test/porousmediumflow/1p/incompressible/main.cc
              LABELS porousmediumflow 1p
              COMPILE_DEFINITIONS TYPETAG=OnePInternalDirichletTpfa
              COMPILE_DEFINITIONS NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_internaldirichlet_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_internaldirichlet_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_internaldirichlet_tpfa ${CMAKE_SOURCE_DIR}/test/porousmediumflow/1p/incompressible/params.input -Problem.Name test_1p_internaldirichlet_tpfa -Problem.EnableGravity false")

dumux_add_test(NAME test_1p_internaldirichlet_box
              SOURCES ${CMAKE_SOURCE_DIR}/test/porousmediumflow/1p/incompressible/main.cc
              LABELS porousmediumflow 1p
              COMPILE_DEFINITIONS TYPETAG=OnePInternalDirichletBox
              COMPILE_DEFINITIONS NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_internaldirichlet_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_internaldirichlet_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_internaldirichlet_box ${CMAKE_SOURCE_DIR}/test/porousmediumflow/1p/incompressible/params.input -Problem.Name test_1p_internaldirichlet_box -Problem.EnableGravity false")
