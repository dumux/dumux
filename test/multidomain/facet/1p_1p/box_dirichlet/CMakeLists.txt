# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "grids" "params.input")

dumux_add_test(NAME test_md_facet_1p1p_box_dirichlet
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1p
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_box_dirichlet_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_box_dirichlet_bulk-00001.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_box_dirichlet_lowdim-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_box_dirichlet_lowdim-00001.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_box_dirichlet params.input")
