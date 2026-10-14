# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input convergencetest.py plot_results.py)

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
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS test_ff_navierstokes_taylorcouette)
