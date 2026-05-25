# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params_mixed.input)
file(CREATE_LINK
    ${CMAKE_SOURCE_DIR}/test/solidmechanics/hyperelastic/cooks_membrane_2d/cooksmembrane.msh
    ${CMAKE_CURRENT_BINARY_DIR}/cooksmembrane.msh SYMBOLIC)

dumux_add_test(
    NAME test_hyperelastic_cooks_membrane_mini
    LABELS geomechanics elastic
    SOURCES main_mixed.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane_mini-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_mini_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_mini params_mixed.input"
)
dumux_add_test(
    NAME test_hyperelastic_cooks_membrane_taylorhood
    LABELS geomechanics elastic
    SOURCES main_taylorhood.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane_taylorhood-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_taylorhood_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_taylorhood params_mixed.input"
)
