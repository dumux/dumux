# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params_mixed3d.input cooksmembrane3d.msh)

dumux_add_test(
    NAME test_hyperelastic_cooks_membrane_mini_3d
    LABELS geomechanics elastic
    SOURCES main_mixed3d.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane_mini3d-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_mini3d_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_mini_3d params_mixed3d.input -SpatialParams.LambdaValues 432.099 -Grid.Refinement 1"
)
dumux_add_test(
    NAME test_hyperelastic_cooks_membrane_taylorhood_3d
    LABELS geomechanics elastic
    SOURCES main_taylorhood3d.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane_taylorhood3d-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_taylorhood3d_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_taylorhood_3d params_mixed3d.input -SpatialParams.LambdaValues 432.099 -Grid.Refinement 1"
)
