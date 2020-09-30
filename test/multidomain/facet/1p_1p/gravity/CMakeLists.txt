dune_symlink_to_source_files(FILES "grids" "params.input")

# test gravity using different values for xi
add_executable(test_md_facet_1p1p_gravity_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_gravity_tpfa
                           PUBLIC BULKTYPETAG=OnePBulkTpfa
                                  LOWDIMTYPETAG=OnePLowDimTpfa
                                  LOWDIMGRIDTYPE=Dune::FoamGrid<1,2>
                                  BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>)

add_executable(test_md_facet_1p1p_gravity_mpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_gravity_mpfa
                           PUBLIC BULKTYPETAG=OnePBulkMpfa
                                  LOWDIMTYPETAG=OnePLowDimMpfa
                                  LOWDIMGRIDTYPE=Dune::FoamGrid<1,2>
                                  BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>)

dumux_add_test(NAME test_md_facet_1p1p_gravity_xi1_tpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_gravity_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_tpfa_bulk-00000.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_tpfa_bulk-00001.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_tpfa_lowdim-00000.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_tpfa_lowdim-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_tpfa params.input \
                                                              -Vtk.OutputName test_md_facet_1p1p_gravity_xi1_tpfa")


dumux_add_test(NAME test_md_facet_1p1p_gravity_xi066_tpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_gravity_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_tpfa_bulk-00000.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_tpfa_bulk-00001.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_tpfa_lowdim-00000.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_tpfa_lowdim-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_tpfa params.input  \
                                                               -FacetCoupling.Xi 0.66 -Vtk.OutputName test_md_facet_1p1p_gravity_xi066_tpfa")



dumux_add_test(NAME test_md_facet_1p1p_gravity_xi1_mpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_gravity_mpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_mpfa_bulk-00000.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_mpfa_bulk-00001.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_mpfa_lowdim-00000.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_mpfa_lowdim-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_mpfa params.input \
                                                              -Vtk.OutputName test_md_facet_1p1p_gravity_xi1_mpfa")

dumux_add_test(NAME test_md_facet_1p1p_gravity_xi066_mpfa
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_gravity_mpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_mpfa_bulk-00000.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_mpfa_bulk-00001.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_mpfa_lowdim-00000.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_mpfa_lowdim-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_mpfa params.input  \
                                                               -FacetCoupling.Xi 0.66
                                                               -Vtk.OutputName test_md_facet_1p1p_gravity_xi066_mpfa")

dumux_add_test(NAME test_md_facet_1p1p_gravity_surface_tpfa
               LABELS multidomain multidomain_facet 1p
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               SOURCES main.cc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkTpfa
               COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePLowDimTpfa
               COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<2,3,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS LOWDIMGRIDTYPE=Dune::FoamGrid<1,3>
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_tpfa_bulk-00000.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_tpfa_bulk-00001.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_tpfa_lowdim-00000.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_tpfa_lowdim-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_tpfa params.input \
                                                               -Vtk.OutputName test_md_facet_1p1p_gravity_surface_tpfa \
                                                               -Grid.File grids/gravity_surface.msh")

dumux_add_test(NAME test_md_facet_1p1p_gravity_surface_mpfa
               LABELS multidomain multidomain_facet 1p
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               SOURCES main.cc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkMpfa
               COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePLowDimMpfa
               COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<2,3,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS LOWDIMGRIDTYPE=Dune::FoamGrid<1,3>
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_mpfa_bulk-00000.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_mpfa_bulk-00001.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_mpfa_lowdim-00000.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_mpfa_lowdim-00001.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_surface_mpfa params.input  \
                                                                -Vtk.OutputName test_md_facet_1p1p_gravity_surface_mpfa \
                                                                -Grid.File grids/gravity_surface.msh")
