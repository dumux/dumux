# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params3d.input cooksmembrane3d.msh)

dumux_add_test(
    NAME test_hyperelastic_cooks_membrane_3d
    LABELS geomechanics elastic
    SOURCES main3d.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane3d-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane3d_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_3d params3d.input"
)
dumux_add_test(
    NAME test_hyperelastic_cooks_membrane_3d_pq2
    LABELS geomechanics elastic
    SOURCES main3d_pq2.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane3d_pq2-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane3d_pq2_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_3d_pq2 params3d.input"
)
