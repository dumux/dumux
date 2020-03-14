dune_symlink_to_source_files(FILES "grids" "params.input")

dumux_add_test(NAME test_md_facet_1pnc1pnc_tpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkTpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetTpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_tpfa_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_tpfa_facet-00018.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_tpfa params.input -Vtk.OutputName test_md_facet_1pnc1pnc_tpfa")

dumux_add_test(NAME test_md_facet_1pnc1pnc_mpfa
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkMpfa
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetTpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_tpfa_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_mpfa_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_tpfa_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_mpfa_facet-00018.vtp
                         --zeroThreshold {"delp":1e-9}
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_mpfa params.input -Vtk.OutputName test_md_facet_1pnc1pnc_mpfa")

dumux_add_test(NAME test_md_facet_1pnc1pnc_box
               SOURCES main.cc
               LABELS multidomain multidomain_facet 1pnc
               COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkBox
               COMPILE_DEFINITIONS FACETTYPETAG=OnePNCFacetBox
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_box_bulk-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box_bulk-00018.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1pnc1pnc_box_facet-reference.vtp
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box_facet-00018.vtp
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1pnc1pnc_box params.input -Vtk.OutputName test_md_facet_1pnc1pnc_box")
