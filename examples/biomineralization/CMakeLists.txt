# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES injection_checkpoints.dat injection_type.dat params.input)

# compile MICP simplified chemistry column setup
dumux_add_test(NAME example_biomineralization
               LABELS porousmediumflow example
               TIMEOUT 3600
               SOURCES main.cc
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_biomineralization-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/example_biomineralization-00018.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/example_biomineralization params.input -Problem.Name  example_biomineralization")
