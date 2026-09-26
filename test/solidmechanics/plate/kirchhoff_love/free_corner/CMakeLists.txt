# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input run_free_corner.py)

dumux_add_test(NAME test_kirchhoff_love_free_corner
              LABELS solidmechanics plate
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND Gmsh_FOUND )"
              COMMAND ./run_free_corner.py
              CMD_ARGS ./test_kirchhoff_love_free_corner)
