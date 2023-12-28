# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# These tests test the behaviour of the box-facet-coupling model in the case that a fracture ends
# at the boundary, or extends beyond the boundary. In both cases, the dofs at the boundary of the
# bulk domain have to be duplicated in order to account for possible pressure jumps. The tests use
# a very low fracture permeability, and we use the same bulk reference solution to test against.

dune_symlink_to_source_files(FILES "grids" "params.input")

add_executable(test_md_facet_1p1p_box_dof_enrichment EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_md_facet_1p1p_box_dof_enrichment_extended_fracture
               TARGET test_md_facet_1p1p_box_dof_enrichment
               LABELS multidomain multidomain_facet 1p
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_dof_enrichment_extended_fracture_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_dof_enrichment_extended_fracture_bulk-00001.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_dof_enrichment_extended_fracture_lowdim-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_dof_enrichment_extended_fracture_lowdim-00001.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_box_dof_enrichment params.input \
                                                                            -Grid.File grids/extended_fracture.msh \
                                                                            -Vtk.OutputName test_md_facet_1p1p_dof_enrichment_extended_fracture")

dumux_add_test(NAME test_md_facet_1p1p_box_dof_enrichment_short_fracture
               TARGET test_md_facet_1p1p_box_dof_enrichment
               LABELS multidomain multidomain_facet 1p
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_dof_enrichment_extended_fracture_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_dof_enrichment_short_fracture_bulk-00001.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_dof_enrichment_short_fracture_lowdim-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_dof_enrichment_short_fracture_lowdim-00001.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_box_dof_enrichment params.input \
                                                                            -Grid.File grids/short_fracture.msh \
                                                                            -Vtk.OutputName test_md_facet_1p1p_dof_enrichment_short_fracture")
