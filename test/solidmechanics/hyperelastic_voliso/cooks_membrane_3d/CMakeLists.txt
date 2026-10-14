# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

foreach(_file cooksmembrane3d.msh params3d.input)
    file(CREATE_LINK
        ${CMAKE_SOURCE_DIR}/test/solidmechanics/hyperelastic/cooks_membrane_3d/${_file}
        ${CMAKE_CURRENT_BINARY_DIR}/${_file} SYMBOLIC)
endforeach()

dumux_add_test(
    NAME test_hyperelastic_cooks_membrane_mini_3d
    LABELS geomechanics elastic
    SOURCES main_mixed3d.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_cooks_membrane_mini3d-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_mini3d_lambda432-00000.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_mini_3d params3d.input"
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
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_cooks_membrane_taylorhood_3d params3d.input"
)
