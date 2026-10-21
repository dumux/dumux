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

# two skewed triangles of a surface grid embedded in 3D
dumux_add_test(NAME test_1p_hydrostatic_surface_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS TYPETAG=OnePHydrostaticSurfaceTpfa
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p_hydrostatic_surface_tpfa
              CMD_ARGS params.input -Grid.File grids/twoskewedtriangles3d.dgf -SpatialParams.PermeabilityJumpX 0.5)

# network with a branching point embedded in 2D
dumux_add_test(NAME test_1p_hydrostatic_network_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS TYPETAG=OnePHydrostaticNetworkTpfa
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p_hydrostatic_network_tpfa
              CMD_ARGS params.input -Grid.File grids/branchingnetwork2d.dgf)
