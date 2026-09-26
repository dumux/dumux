# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_submanager.input)

dumux_add_test(NAME test_kirchhoff_love_plate
              THREADS 4
              LABELS solidmechanics plate
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND Gmsh_FOUND )"
              COMMAND ../convergencetest.py
              CMD_ARGS test_kirchhoff_love_plate)

dumux_add_test(NAME test_kirchhoff_love_plate_submanager
              LABELS solidmechanics plate
              SOURCES test_submanager.cc
              CMAKE_GUARD dune-foamgrid_FOUND)

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
