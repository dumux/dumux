# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_md_embedded_1d3d_lowdimpartition
               LABELS multidomain multidomain_embedded parallel
               SOURCES main.cc
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 300)
