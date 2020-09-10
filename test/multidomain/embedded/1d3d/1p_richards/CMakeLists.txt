dumux_add_test(NAME test_md_embedded_1d3d_1p_richards_tpfatpfa
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfatpfa_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_tpfatpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfatpfa_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfatpfa params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p_richards_tpfatpfa")

dumux_add_test(NAME test_md_embedded_1d3d_1p_richards_tpfabox
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfabox_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_tpfabox_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfabox_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfabox params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p_richards_tpfabox")
dune_symlink_to_source_files(FILES "params.input")
