dune_symlink_to_source_files(FILES "grids" "params.input")

dumux_add_test(NAME test_md_facet_1pnc1pnc_tpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkTpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetTpfa
               COMPILE_DEFINITIONS DIMWORLD=2
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_tpfa_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_tpfa_facet-00018.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_tpfa params.input -Vtk.OutputName test_md_facet_1pnc1pnc_tpfa")

dumux_add_test(NAME test_md_facet_1pnc1pnc_mpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkMpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetTpfa
               COMPILE_DEFINITIONS DIMWORLD=2
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_mpfa_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_mpfa_facet-00018.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_mpfa params.input -Vtk.OutputName test_md_facet_1pnc1pnc_mpfa")

dumux_add_test(NAME test_md_facet_1pnc1pnc_box
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkBox
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetBox
               COMPILE_DEFINITIONS DIMWORLD=2
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_box_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_box_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box_facet-00018.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box params.input -Vtk.OutputName test_md_facet_1pnc1pnc_box")

dumux_add_test(NAME test_md_facet_1pnc1pnc_box_mixedbcs
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkBox
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetBox
               COMPILE_DEFINITIONS DIMWORLD=2
               COMPILE_DEFINITIONS USEMIXEDBCS=true
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_box_mixedbcs_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box_mixedbcs_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_box_mixedbcs_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box_mixedbcs_facet-00018.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box_mixedbcs params.input -Vtk.OutputName test_md_facet_1pnc1pnc_box_mixedbcs")

dumux_add_test(NAME test_md_facet_1pncni1pncni_tpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCNIBulkTpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCNIFacetTpfa
               COMPILE_DEFINITIONS DIMWORLD=2
               COMPILE_DEFINITIONS BULKGRID=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_tpfa_bulk-00019.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_tpfa_facet-00019.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_tpfa params.input -Vtk.OutputName test_md_facet_1pncni1pncni_tpfa")

dumux_add_test(NAME test_md_facet_1pncni1pncni_mpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCNIBulkMpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCNIFacetTpfa
               COMPILE_DEFINITIONS DIMWORLD=2
               COMPILE_DEFINITIONS BULKGRID=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_mpfa_bulk-00019.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_mpfa_facet-00019.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_mpfa params.input -Vtk.OutputName test_md_facet_1pncni1pncni_mpfa")

dumux_add_test(NAME test_md_facet_1pncni1pncni_box
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCNIBulkBox
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCNIFacetBox
               COMPILE_DEFINITIONS DIMWORLD=2
               COMPILE_DEFINITIONS BULKGRID=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_box_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_box_bulk-00019.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_box_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_box_facet-00019.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_box params.input -Vtk.OutputName test_md_facet_1pncni1pncni_box")


dumux_add_test(NAME test_md_facet_1pnc1pnc_surface_tpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkTpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetTpfa
               COMPILE_DEFINITIONS DIMWORLD=3
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_surface_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_tpfa_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_surface_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_tpfa_facet-00018.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_tpfa params.input -Vtk.OutputName test_md_facet_1pnc1pnc_surface_tpfa -Grid.File grids/square_surface.msh")

dumux_add_test(NAME test_md_facet_1pnc1pnc_surface_mpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkMpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetTpfa
               COMPILE_DEFINITIONS DIMWORLD=3
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --zeroThreshold {"delp":1e-8}
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_surface_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_mpfa_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_surface_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_mpfa_facet-00018.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_mpfa params.input -Vtk.OutputName test_md_facet_1pnc1pnc_surface_mpfa -Grid.File grids/square_surface.msh")

dumux_add_test(NAME test_md_facet_1pnc1pnc_surface_box
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkBox
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetBox
               COMPILE_DEFINITIONS DIMWORLD=3
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_surface_box_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_box_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_surface_box_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_box_facet-00018.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_surface_box params.input -Vtk.OutputName test_md_facet_1pnc1pnc_surface_box -Grid.File grids/square_surface.msh")

dumux_add_test(NAME test_md_facet_1pncni1pncni_surface_tpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCNIBulkTpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCNIFacetTpfa
               COMPILE_DEFINITIONS DIMWORLD=3
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_surface_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_tpfa_bulk-00019.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_surface_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_tpfa_facet-00019.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_tpfa params.input -Vtk.OutputName test_md_facet_1pncni1pncni_surface_tpfa -Grid.File grids/square_surface.msh")

dumux_add_test(NAME test_md_facet_1pncni1pncni_surface_mpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCNIBulkMpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCNIFacetTpfa
               COMPILE_DEFINITIONS DIMWORLD=3
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --zeroThreshold {"delp":1e-8}
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_surface_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_mpfa_bulk-00020.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_surface_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_mpfa_facet-00020.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_mpfa params.input -Vtk.OutputName test_md_facet_1pncni1pncni_surface_mpfa -Grid.File grids/square_surface.msh")

dumux_add_test(NAME test_md_facet_1pncni1pncni_surface_box
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCNIBulkBox
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCNIFacetBox
               COMPILE_DEFINITIONS DIMWORLD=3
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_surface_box_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_box_bulk-00019.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pncni1pncni_surface_box_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_box_facet-00019.vtp
                         --zeroThreshold {"delp":1e-3}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pncni1pncni_surface_box params.input -Vtk.OutputName test_md_facet_1pncni1pncni_surface_box -Grid.File grids/square_surface.msh")
