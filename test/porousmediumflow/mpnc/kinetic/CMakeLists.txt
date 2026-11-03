# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# build target for the full kinetic test problem
dumux_add_test(NAME test_mpnc_kinetic_box
              LABELS porousmediumflow mpnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=EvaporationAtmosphereBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_kinetic_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_kinetic_box-00011.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_kinetic_box params.input -Problem.Name test_mpnc_kinetic_box")
