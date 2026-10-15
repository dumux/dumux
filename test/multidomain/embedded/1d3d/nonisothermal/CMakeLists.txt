# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_wellbore_heat_transport_surface
              SOURCES main.cc
              COMPILE_DEFINITIONS COUPLINGMODE=Embedded1d3dCouplingMode::Surface
              COMMAND ${DUMUX_RUNTEST}
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_surface_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_surface_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_surface_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_surface_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_surface -TimeLoop.EpisodeLength 432000 \
                                   -Vtk.OutputName test_wellbore_heat_transport_surface -MixedDimension.NumCircleSegments 8")

dumux_add_test(NAME test_wellbore_heat_transport_kernel
              SOURCES main.cc
              COMPILE_DEFINITIONS COUPLINGMODE=Embedded1d3dCouplingMode::Kernel COORDINATETYPE=EquidistantOffsetCoordinates
              COMMAND ${DUMUX_RUNTEST}
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_kernel_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_kernel_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_wellbore_heat_transport_kernel_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_kernel_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_wellbore_heat_transport_kernel -TimeLoop.EpisodeLength 432000 \
                                   -Vtk.OutputName test_wellbore_heat_transport_kernel")

add_input_file_links()
dune_symlink_to_source_files(FILES "grids")
