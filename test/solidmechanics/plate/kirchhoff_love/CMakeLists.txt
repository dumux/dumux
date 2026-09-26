# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_kirchhoff_love_plate
              LABELS solidmechanics plate
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND Gmsh_FOUND )"
              COMMAND ../convergencetest.py
              CMD_ARGS test_kirchhoff_love_plate)

dumux_add_test(NAME test_kirchhoff_love_plate_quadrature
              LABELS solidmechanics plate
              SOURCES main.cc
              COMPILE_DEFINITIONS SCVF_QUADRATURE_ORDER=3
              CMAKE_GUARD "( dune-foamgrid_FOUND AND Gmsh_FOUND )"
              COMMAND ../convergencetest.py
              CMD_ARGS test_kirchhoff_love_plate_quadrature)

set_tests_properties(test_kirchhoff_love_plate test_kirchhoff_love_plate_quadrature
                     PROPERTIES RESOURCE_LOCK kirchhoff_love_plate_mesh)

add_subdirectory(annulus)
add_subdirectory(annulus_clamped)
add_subdirectory(free_corner)
add_subdirectory(navier_square)
