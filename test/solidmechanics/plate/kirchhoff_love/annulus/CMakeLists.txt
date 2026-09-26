# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_kirchhoff_love_plate_annulus
              LABELS solidmechanics plate
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND Gmsh_FOUND )"
              COMMAND ../../convergencetest.py
              CMD_ARGS test_kirchhoff_love_plate_annulus --geo ../../annulus.geo --mesh annulus.msh
                       --clmax 0.125 0.0625 0.03125 0.015625)
