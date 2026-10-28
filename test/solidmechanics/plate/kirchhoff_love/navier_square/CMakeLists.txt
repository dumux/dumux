# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input run_navier_square.py run_clamped_square.py)

dumux_add_test(NAME test_kirchhoff_love_navier_square
              LABELS solidmechanics plate
              SOURCES main.cc
              CMAKE_GUARD "( dune-uggrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ./run_navier_square.py
              CMD_ARGS ./test_kirchhoff_love_navier_square)

# all four edges clamped, where the curl potential is not constant, against a Legendre-Galerkin reference
dumux_add_test(NAME test_kirchhoff_love_clamped_square
              LABELS solidmechanics plate
              TARGET test_kirchhoff_love_navier_square
              CMAKE_GUARD "( dune-uggrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ./run_clamped_square.py
              CMD_ARGS ./test_kirchhoff_love_navier_square)
