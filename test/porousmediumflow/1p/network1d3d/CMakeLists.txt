# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dune_symlink_to_source_files(FILES grids)

# dim < dimWorld tests with Dune::Foamgrid<1,3>
dumux_add_test(NAME test_1p_network1d3d_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TubesTestCCTpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_network1d3d_tpfa params.input -Problem.Name test_1p_network1d3d_tpfa"
                       --levels 0 1 2 3 4 5
                       --level-args "-Grid.Refinement {}"
                       --errors test_1p_network1d3d_tpfa-errors.csv
                       --check "L2(p) all 1.5 2.5")

dumux_add_test(NAME test_1p_network1d3d_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TubesTestBox
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_network1d3d_box params.input -Problem.Name test_1p_network1d3d_box"
                       --levels 0 1 2 3 4 5
                       --level-args "-Grid.Refinement {}"
                       --errors test_1p_network1d3d_box-errors.csv
                       --check "L2(p) all 1.5 2.5")
