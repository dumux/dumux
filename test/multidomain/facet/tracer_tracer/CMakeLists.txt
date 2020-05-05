dune_symlink_to_source_files(FILES "grids" "params.input")

dumux_add_test(NAME test_md_facet_tracertracer_tpfa
              LABELS multidomain multidomain_facet tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS ONEPBULKTYPETAG=OnePBulkTpfa
              COMPILE_DEFINITIONS ONEPLOWDIMTYPETAG=OnePLowDimTpfa
              COMPILE_DEFINITIONS TRACERBULKTYPETAG=TracerBulkTpfa
              COMPILE_DEFINITIONS TRACERLOWDIMTYPETAG=TracerLowDimTpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_tpfa_onep_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_tpfa_onep_bulk-00000.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_tpfa_onep_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_tpfa_onep_lowdim-00000.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_tpfa_tracer_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_tpfa_tracer_bulk-00010.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_tpfa_tracer_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_tpfa_tracer_lowdim-00010.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_tpfa params.input \
                                                              -Vtk.OutputName test_md_facet_tracertracer_tpfa \
                                                              -TimeLoop.TEnd 5e4")

dumux_add_test(NAME test_md_facet_tracertracer_mpfa
              LABELS multidomain multidomain_facet tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS ONEPBULKTYPETAG=OnePBulkMpfa
              COMPILE_DEFINITIONS ONEPLOWDIMTYPETAG=OnePLowDimMpfa
              COMPILE_DEFINITIONS TRACERBULKTYPETAG=TracerBulkMpfa
              COMPILE_DEFINITIONS TRACERLOWDIMTYPETAG=TracerLowDimMpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TIMEOUT 1500
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_mpfa_onep_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_mpfa_onep_bulk-00000.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_mpfa_onep_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_mpfa_onep_lowdim-00000.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_mpfa_tracer_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_mpfa_tracer_bulk-00010.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_mpfa_tracer_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_mpfa_tracer_lowdim-00010.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_mpfa params.input \
                                                              -Vtk.OutputName test_md_facet_tracertracer_mpfa \
                                                              -TimeLoop.TEnd 5e4")

dumux_add_test(NAME test_md_facet_tracertracer_box
              LABELS multidomain multidomain_facet tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS ONEPBULKTYPETAG=OnePBulkBox
              COMPILE_DEFINITIONS ONEPLOWDIMTYPETAG=OnePLowDimBox
              COMPILE_DEFINITIONS TRACERBULKTYPETAG=TracerBulkBox
              COMPILE_DEFINITIONS TRACERLOWDIMTYPETAG=TracerLowDimBox
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_box_onep_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_box_onep_bulk-00000.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_box_onep_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_box_onep_lowdim-00000.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_box_tracer_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_box_tracer_bulk-00010.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_tracertracer_box_tracer_lowdim-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_box_tracer_lowdim-00010.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_tracertracer_box params.input \
                                                              -Vtk.OutputName test_md_facet_tracertracer_box \
                                                              -TimeLoop.TEnd 5e4")
