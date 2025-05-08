# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "grids" "params.input")

dumux_add_test(NAME test_md_facet_1p1p_boundary_coupling_tpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1p
               COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkTpfa
               COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePLowDimTpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_boundary_coupling_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_tpfa_bulk-00001.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_boundary_coupling_tpfa_lowdim-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_tpfa_lowdim-00001.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_tpfa params.input \
                                                            -Vtk.OutputName test_md_facet_1p1p_boundary_coupling_tpfa")

dumux_add_test(NAME test_md_facet_1p1p_boundary_coupling_mpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1p
               COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkMpfa
               COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePLowDimMpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_boundary_coupling_mpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_mpfa_bulk-00001.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_boundary_coupling_mpfa_lowdim-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_mpfa_lowdim-00001.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_mpfa params.input \
                                                            -Vtk.OutputName test_md_facet_1p1p_boundary_coupling_mpfa")

dumux_add_test(NAME test_md_facet_1p1p_boundary_coupling_box
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1p
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkBox
               COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePLowDimBox
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_boundary_coupling_box_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_box_bulk-00001.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_boundary_coupling_box_lowdim-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_box_lowdim-00001.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_box params.input \
                                                            -Vtk.OutputName test_md_facet_1p1p_boundary_coupling_box")

dumux_add_test(NAME test_md_facet_1p1p_boundary_coupling_box_caching
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1p
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkBox
               COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePLowDimBox
               COMPILE_DEFINITIONS USEBULKGGCACHING=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_boundary_coupling_box_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_box_caching_bulk-00001.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_boundary_coupling_box_lowdim-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_box_caching_lowdim-00001.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_boundary_coupling_box_caching params.input \
                                                            -Vtk.OutputName test_md_facet_1p1p_boundary_coupling_box_caching")
