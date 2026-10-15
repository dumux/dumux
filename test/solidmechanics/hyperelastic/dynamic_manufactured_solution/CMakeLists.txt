# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input convergencetest.py)

dumux_add_test(NAME test_hyperelastic_dynamic_manufactured_solution_pq2fe
               LABELS solidmechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=HyperelasticDynamicMMSPQ2FE
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ./convergencetest.py
               CMD_ARGS test_hyperelastic_dynamic_manufactured_solution_pq2fe 2.8 1.8)

dumux_add_test(NAME test_hyperelastic_dynamic_manufactured_solution_pq2hybrid
               LABELS solidmechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=HyperelasticDynamicMMSPQ2Hybrid
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ./convergencetest.py
               CMD_ARGS test_hyperelastic_dynamic_manufactured_solution_pq2hybrid 1.8 1.8)
