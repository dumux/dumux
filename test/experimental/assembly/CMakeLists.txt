# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_experimental_cvfe_analytic_jacobian
               LABELS unit assembly
               SOURCES main.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_experimental_cvfe_analytic_jacobian)

dumux_add_test(NAME test_experimental_assembler_tpfa
               LABELS unit assembly
               SOURCES tpfa.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_experimental_assembler_tpfa)

dumux_add_test(NAME test_experimental_assembler_fcstaggered
               LABELS unit assembly
               SOURCES fcstaggered.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_experimental_assembler_fcstaggered)
