dune_symlink_to_source_files(FILES "grids" "params.input")

add_executable(test_md_facet_1p1p_linearprofile_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_linearprofile_tpfa
                           PUBLIC BULKTYPETAG=OnePBulkTpfa
                                  LOWDIMTYPETAG=OnePLowDimTpfa
                                  LOWDIMGRIDTYPE=Dune::FoamGrid<1,2>
                                  BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>)

add_executable(test_md_facet_1p1p_linearprofile_surface_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_linearprofile_surface_tpfa
                           PUBLIC BULKTYPETAG=OnePBulkTpfa
                                  LOWDIMTYPETAG=OnePLowDimTpfa
                                  LOWDIMGRIDTYPE=Dune::FoamGrid<1,3>
                                  BULKGRIDTYPE=Dune::ALUGrid<2,3,Dune::cube,Dune::nonconforming>)

add_executable(test_md_facet_1p1p_linearprofile_mpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_linearprofile_mpfa
                           PUBLIC BULKTYPETAG=OnePBulkMpfa
                                  LOWDIMTYPETAG=OnePLowDimMpfa
                                  LOWDIMGRIDTYPE=Dune::FoamGrid<1,2>
                                  BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>)

add_executable(test_md_facet_1p1p_linearprofile_surface_mpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_linearprofile_surface_mpfa
                            PUBLIC BULKTYPETAG=OnePBulkMpfa
                                   LOWDIMTYPETAG=OnePLowDimMpfa
                                   LOWDIMGRIDTYPE=Dune::FoamGrid<1,3>
                                   BULKGRIDTYPE=Dune::ALUGrid<2,3,Dune::cube,Dune::nonconforming>)

dumux_add_test(NAME test_md_facet_1p1p_linearprofile_xi1_tpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_linearprofile_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_tpfa_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_xi1_tpfa_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_tpfa_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_xi1_tpfa_lowdim-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_tpfa params.input")

dumux_add_test(NAME test_md_facet_1p1p_linearprofile_xi066_tpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_linearprofile_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_tpfa_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_xi066_tpfa_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_tpfa_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_xi066_tpfa_lowdim-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_tpfa params.input \
                                                               -FacetCoupling.Xi 0.66 \
                                                               -Vtk.OutputName test_md_facet_1p1p_linearprofile_xi066_tpfa")

dumux_add_test(NAME test_md_facet_1p1p_linearprofile_surface_xi1_tpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_linearprofile_surface_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_surface_tpfa_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_surface_xi1_tpfa_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_surface_tpfa_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_surface_xi1_tpfa_lowdim-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_surface_tpfa params.input \
                                                              -Grid.File grids/linear_surface.msh \
                                                              -Vtk.OutputName test_md_facet_1p1p_linearprofile_surface_xi1_tpfa")

dumux_add_test(NAME test_md_facet_1p1p_linearprofile_xi1_mpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_linearprofile_mpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_tpfa_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_xi1_mpfa_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_tpfa_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_xi1_mpfa_lowdim-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_mpfa params.input \
                                                              -Vtk.OutputName test_md_facet_1p1p_linearprofile_xi1_mpfa")

dumux_add_test(NAME test_md_facet_1p1p_linearprofile_xi066_mpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_linearprofile_mpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_tpfa_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_xi066_mpfa_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_tpfa_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_xi066_mpfa_lowdim-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_mpfa params.input \
                                                            -FacetCoupling.Xi 0.66 \
                                                            -Vtk.OutputName test_md_facet_1p1p_linearprofile_xi066_mpfa")

dumux_add_test(NAME test_md_facet_1p1p_linearprofile_surface_xi1_mpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_linearprofile_surface_mpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_surface_tpfa_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_surface_xi1_mpfa_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_linearprofile_surface_tpfa_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_surface_xi1_mpfa_lowdim-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_linearprofile_surface_mpfa params.input \
                                                              -Grid.File grids/linear_surface.msh \
                                                              -Vtk.OutputName test_md_facet_1p1p_linearprofile_surface_xi1_mpfa")
