# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" grids)

# hydrostatic pressure on two skewed cells (non K-orthogonal), self-checking
dumux_add_test(NAME test_1p_hydrostatic_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS TYPETAG=OnePHydrostaticTpfa
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p_hydrostatic_tpfa
              CMD_ARGS params.input)

dumux_add_test(NAME test_1p_hydrostatic_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS TYPETAG=OnePHydrostaticBox
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p_hydrostatic_box
              CMD_ARGS params.input)
