# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input cooksmembrane.msh)

dumux_add_test(
    NAME test_cooks_membrane_box
    LABELS geomechanics elastic
    SOURCES main.cc
    COMPILE_DEFINITIONS DISCRETIZATION_MODEL=BoxModel
    CMAKE_GUARD "dune-alugrid_FOUND"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_cooks_membrane_box-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_cooks_membrane_box-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_cooks_membrane_box params.input"
)

dumux_add_test(
    NAME test_cooks_membrane_pq1bubble
    LABELS geomechanics elastic
    SOURCES main.cc
    COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ1BubbleModel
    CMAKE_GUARD "dune-alugrid_FOUND"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_cooks_membrane_pq1bubble-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_cooks_membrane_pq1bubble-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_cooks_membrane_pq1bubble params.input -Problem.Name test_cooks_membrane_pq1bubble"
)

dumux_add_test(
    NAME test_cooks_membrane_pq2
    LABELS geomechanics elastic
    SOURCES main.cc
    COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ2HybridModel
    CMAKE_GUARD "dune-alugrid_FOUND"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_cooks_membrane_pq2-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_cooks_membrane_pq2-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_cooks_membrane_pq2 params.input -Problem.Name test_cooks_membrane_pq2"
)
