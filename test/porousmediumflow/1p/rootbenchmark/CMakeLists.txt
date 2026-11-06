# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_1p_rootbenchmark_tpfa
               LABELS porousmediumflow 1p
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=RootBenchmarkCCTpfa
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_rootbenchmark_tpfa params.input -Problem.EnableGravity false -Problem.Name test_1p_rootbenchmark_tpfa"
                        --levels 0 1 2 3 4 5
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_rootbenchmark_tpfa-errors.csv
                        --check "L2Rel(p) all 1.5 2.5")

dumux_add_test(NAME test_1p_rootbenchmark_tpfa_gravity
               TARGET test_1p_rootbenchmark_tpfa
               LABELS porousmediumflow 1p
               COMPILE_DEFINITIONS TYPETAG=RootBenchmarkCCTpfa
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_rootbenchmark_tpfa params.input -Problem.Name test_1p_rootbenchmark_tpfa_gravity"
                        --levels 0 1 2 3 4 5
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_rootbenchmark_tpfa_gravity-errors.csv
                        --check "L2Rel(p) all 1.5 2.5")
