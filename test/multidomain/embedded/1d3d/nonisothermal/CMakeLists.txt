# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_wellbore_heat_transport_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=SoilCC
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/dumux_runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_tpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_tpfa_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_tpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_tpfa_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_tpfa -TimeLoop.TEnd 7200 \
                                   -MixedDimension.NumCircleSegments 8 -Vtk.OutputName test_wellbore_heat_transport_tpfa")

dumux_add_test(NAME test_wellbore_heat_transport_box
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=SoilBox
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/dumux_runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_box_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_box_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_box_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_box_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_box -TimeLoop.TEnd 7200 \
                                   -MixedDimension.NumCircleSegments 8 -Vtk.OutputName test_wellbore_heat_transport_box")


# The kernel (distributed source) method smears the exchange flux of a wellbore segment over
# the bulk elements overlapping the kernel support, which requires an equidistant bulk grid.
dumux_add_test(NAME test_wellbore_heat_transport_kernel
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=SoilCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel COORDINATETYPE=EquidistantOffsetCoordinates
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/dumux_runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_kernel_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_kernel_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_kernel_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_kernel_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_kernel params_kernel.input -TimeLoop.TEnd 7200 \
                                   -Soil.Grid.Cells \"20 20 10\" -Vtk.OutputName test_wellbore_heat_transport_kernel")

add_input_file_links()
dune_symlink_to_source_files(FILES "grids" "convergence_test.py" "convergence_test_kernel.py" "ramey.py" "ogs_reference.py")
