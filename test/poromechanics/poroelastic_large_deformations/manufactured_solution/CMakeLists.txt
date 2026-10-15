# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input convergencetest.py)

dumux_add_test(NAME test_poroelastic_large_def_manufactured_solution
               LABELS poromechanics
               SOURCES main.cc
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ./convergencetest.py
               CMD_ARGS test_poroelastic_large_def_manufactured_solution)

dumux_add_test(NAME test_poroelastic_large_def_manufactured_solution_pq2fe
               LABELS poromechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS MOMENTUM_PQ2FE=1
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ./convergencetest.py
               CMD_ARGS test_poroelastic_large_def_manufactured_solution_pq2fe)
