# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_parallel_indexset_empty_process
               SOURCES test_parallel_indexset_empty_process.cc
               CMAKE_GUARD dune-uggrid_FOUND
               MPI_RANKS 1 2 4
               TIMEOUT 1200)
