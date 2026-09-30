# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_submanager.input)

dumux_add_test(NAME test_kirchhoff_love_plate
              LABELS solidmechanics plate
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND Gmsh_FOUND )"
              COMMAND ../convergencetest.py
              CMD_ARGS test_kirchhoff_love_plate)

dumux_add_test(NAME test_kirchhoff_love_plate_submanager
              LABELS solidmechanics plate
              SOURCES test_submanager.cc
              CMAKE_GUARD dune-foamgrid_FOUND)
