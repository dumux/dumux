# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES ball.msh rectangle.msh)
dumux_add_test(NAME test_parallel_partition_scotch
               SOURCES main.cc
               CMAKE_GUARD "( dune-uggrid_FOUND AND PTSCOTCH_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200)
