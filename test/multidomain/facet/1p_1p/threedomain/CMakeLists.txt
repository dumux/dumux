dune_symlink_to_source_files(FILES "grids" "params.input")

dumux_add_test(NAME test_md_facet_1p1p_threedomain_tpfa
              LABELS multidomain multidomain_facet 1p
              COMPILE_DEFINITIONS BULKTYPETAG=Dumux::Properties::TTag::OnePBulkTpfa
              COMPILE_DEFINITIONS FACETTYPETAG=Dumux::Properties::TTag::OnePFacetTpfa
              COMPILE_DEFINITIONS EDGETYPETAG=Dumux::Properties::TTag::OnePEdgeTpfa
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_tpfa_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_tpfa_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_tpfa_facet-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_tpfa_facet-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_tpfa_edge-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_tpfa_edge-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_tpfa params.input
                       -Vtk.OutputName test_md_facet_1p1p_threedomain_tpfa")

# DO not test mpfa yet because it is too inefficient.
# As soon as grid-wide flux variables caching is supported in multidomain
# simulations we can probably include this test as it will speed it up significantly
# TODO create reference solutions
# dumux_add_test(NAME test_md_facet_1p1p_threedomain_mpfa
#               LABELS multidomain multidomain_facet 1p
#               COMPILE_DEFINITIONS BULKTYPETAG=Dumux::Properties::TTag::OnePBulkMpfa
#               COMPILE_DEFINITIONS FACETTYPETAG=Dumux::Properties::TTag::OnePFacetMpfa
#               COMPILE_DEFINITIONS EDGETYPETAG=Dumux::Properties::TTag::OnePEdgeMpfa
#               SOURCES main.cc
#               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
#               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#               CMD_ARGS  --script fuzzy
#                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_mpfa_bulk-reference.vtu
#                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_mpfa_bulk-00001.vtu
#                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_tpfa_facet-reference.vtu
#                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_mpfa_facet-00001.vtu
#                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_mpfa_edge-reference.vtp
#                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_mpfa_edge-00001.vtp
#                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_mpfa params.input
#                        -Vtk.OutputName test_md_facet_1p1p_threedomain_mpfa")

dumux_add_test(NAME test_md_facet_1p1p_threedomain_box
              LABELS multidomain multidomain_facet 1p
              COMPILE_DEFINITIONS BULKTYPETAG=Dumux::Properties::TTag::OnePBulkBox
              COMPILE_DEFINITIONS FACETTYPETAG=Dumux::Properties::TTag::OnePFacetBox
              COMPILE_DEFINITIONS EDGETYPETAG=Dumux::Properties::TTag::OnePEdgeBox
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_box_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_box_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_box_facet-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_box_facet-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_box_edge-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_box_edge-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_box params.input
                       -Vtk.OutputName test_md_facet_1p1p_threedomain_box")
