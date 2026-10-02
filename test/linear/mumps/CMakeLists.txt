# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_mumps_sequential
               SOURCES main.cc
               CMAKE_GUARD MUMPS_FOUND
               LABELS unit linear)

dumux_add_test(NAME test_mumps_gridadaption
               SOURCES test_mumps_gridadaption.cc
               CMAKE_GUARD MUMPS_FOUND
               MPI_RANKS 1 2 4
               TIMEOUT 300
               LABELS unit linear parallel)
