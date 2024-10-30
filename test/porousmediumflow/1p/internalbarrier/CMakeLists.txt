# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


dumux_add_test(NAME test_1p_internalbarrier_tpfa
              SOURCES ${CMAKE_SOURCE_DIR}/test/porousmediumflow/1p/incompressible/main.cc
              LABELS porousmediumflow 1p
              COMPILE_DEFINITIONS TYPETAG=OnePInternalBarrierTpfa
              COMPILE_DEFINITIONS NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_internalbarrier_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_internalbarrier_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_internalbarrier_tpfa ${CMAKE_SOURCE_DIR}/test/porousmediumflow/1p/incompressible/params.input -Problem.Name test_1p_internaldirichlet_tpfa -Problem.EnableGravity false")
