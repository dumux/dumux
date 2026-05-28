# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(
        NAME test_hyperelastic_incompressible_block_p1b_p1
    LABELS geomechanics elastic
    SOURCES main.cc
    COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION=IncompressibleBlockMomentumP1B
    COMPILE_DEFINITIONS PRESSURE_DISCRETIZATION=IncompressibleBlockPressure
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_incompressible_block_p1b_p1-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_incompressible_block_p1b_p1-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_incompressible_block_p1b_p1 params.input -Problem.Name test_hyperelastic_incompressible_block_p1b_p1"
)

dumux_add_test(
    NAME test_hyperelastic_incompressible_block_p2_p1
    LABELS geomechanics elastic
    SOURCES main.cc
    COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION=IncompressibleBlockMomentumPQ2
    COMPILE_DEFINITIONS PRESSURE_DISCRETIZATION=IncompressibleBlockPressure
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_incompressible_block_p2_p1-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_incompressible_block_p2_p1-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_incompressible_block_p2_p1 params.input -Problem.Name test_hyperelastic_incompressible_block_p2_p1"
)

dumux_add_test(
    NAME test_hyperelastic_incompressible_block_q2_q1
    LABELS geomechanics elastic
    SOURCES main.cc
    COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION=IncompressibleBlockMomentumPQ2
    COMPILE_DEFINITIONS PRESSURE_DISCRETIZATION=IncompressibleBlockPressure
    COMPILE_DEFINITIONS ALU_ELEMENT_TYPE=Dune::cube
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_incompressible_block_q2_q1-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_incompressible_block_q2_q1-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_incompressible_block_q2_q1 params.input -Problem.Name test_hyperelastic_incompressible_block_q2_q1"
)
