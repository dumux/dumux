# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_tracer_multiphase_tpfa
              LABELS porousmediumflow tracer
              TIMEOUT 3600
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_multiphase_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_multiphase_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_multiphase_tpfa params.input -Problem.Name test_tracer_multiphase_tpfa")

dumux_add_test(NAME test_tracer_multiphase_mpfa
              LABELS porousmediumflow tracer
              TIMEOUT 3600
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestMpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_multiphase_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_multiphase_mpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_multiphase_mpfa params.input -Problem.Name test_tracer_multiphase_mpfa")
