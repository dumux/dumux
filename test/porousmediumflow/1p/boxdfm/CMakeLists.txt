# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" "grids")

dumux_add_test(NAME test_1p_boxdfm
              LABELS porousmediumflow 1p boxdfm
              SOURCES main.cc
              CMAKE_GUARD dune-alugrid_FOUND
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_boxdfm-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_boxdfm-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_boxdfm params.input -Problem.Name test_1p_boxdfm -Grid.File grids/durlofsky.msh")
