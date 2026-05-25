# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input incompressible_block.dgf)
dumux_add_test(
    NAME test_hyperelastic_incompressible_block
    LABELS geomechanics elastic
    SOURCES main.cc
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_hyperelastic_incompressible_block-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_incompressible_block-00010.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_incompressible_block params.input -LoadStepping.NSteps 10"
)
