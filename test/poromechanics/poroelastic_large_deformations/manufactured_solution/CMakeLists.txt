# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_poroelastic_large_def_manufactured_solution
               LABELS poromechanics
               SOURCES main.cc
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_large_def_manufactured_solution params.input -Problem.Name test_poroelastic_large_def_manufactured_solution"
                        --levels 8 16 32
                        --level-args "-Grid.Cells '{0} {0}'"
                        --errors test_poroelastic_large_def_manufactured_solution-errors.csv
                        --check "L2(u) all 1.8"
                        --check "L2(p_s) all 1.8"
                        --check "L2(p_f) all 1.8")

dumux_add_test(NAME test_poroelastic_large_def_manufactured_solution_pq2fe
               LABELS poromechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS MOMENTUM_PQ2FE=1
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_large_def_manufactured_solution_pq2fe params.input -Problem.Name test_poroelastic_large_def_manufactured_solution_pq2fe"
                        --levels 8 16 32
                        --level-args "-Grid.Cells '{0} {0}'"
                        --errors test_poroelastic_large_def_manufactured_solution_pq2fe-errors.csv
                        --check "L2(u) all 1.8"
                        --check "L2(p_s) all 1.8"
                        --check "L2(p_f) all 1.8")
