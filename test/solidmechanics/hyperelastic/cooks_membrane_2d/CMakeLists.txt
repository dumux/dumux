# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input cooksmembrane.msh)
file(CREATE_LINK
    ${CMAKE_SOURCE_DIR}/test/solidmechanics/hyperelastic_voliso/cooks_membrane_2d/params_mixed.input
    ${CMAKE_CURRENT_BINARY_DIR}/params_mixed.input SYMBOLIC)

dumux_add_test(
    NAME test_hyperelastic_cooks_membrane
    LABELS geomechanics elastic
    SOURCES main.cc
    CMAKE_GUARD "dune-alugrid_FOUND"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane params.input"
)
dumux_add_test(
    NAME test_hyperelastic_cooks_membrane_pq2
    LABELS geomechanics elastic
    SOURCES main_pq2.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane_pq2-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_pq2_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_pq2 params_mixed.input"
)
