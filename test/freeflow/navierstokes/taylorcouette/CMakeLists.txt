# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input plot_results.py)

dumux_add_test(NAME test_ff_navierstokes_taylorcouette
               LABELS freeflow navierstokes
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_taylorcouette-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_taylorcouette-00000.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_taylorcouette params.input")

dumux_add_test(NAME test_ff_navierstokes_taylorcouette_convergence
               LABELS freeflow navierstokes
               TARGET test_ff_navierstokes_taylorcouette
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_taylorcouette params.input -Problem.Name test_ff_navierstokes_taylorcouette_convergence"
                        --levels 0 1 2
                        --level-args "-Grid.Refinement {}"
                        --errors test_ff_navierstokes_taylorcouette_convergence-errors.csv
                        --check "L2Abs(p) mean 1.9"
                        --check "L2Abs(v) mean 1.9")
