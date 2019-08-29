dune_add_test(NAME test_md_embedded1d3d_1p_richards_tpfatpfa
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilCC
              COMPILE_DEFINITIONS DIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfatpfa_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfatpfa params.input \
                                   -Vtk.OutputName test_md_embedded1d3d_1p_richards_tpfatpfa")

dune_add_test(NAME test_md_embedded1d3d_1p_richards_tpfatpfa_anadiff
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilCC
              COMPILE_DEFINITIONS DIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfatpfa_anadiff_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfatpfa_anadiff_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfatpfa_anadiff params.input \
                                   -Vtk.OutputName test_md_embedded1d3d_1p_richards_tpfatpfa_anadiff")

dune_add_test(NAME test_md_embedded1d3d_1p_richards_tpfabox
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilBox
              COMPILE_DEFINITIONS DIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfabox_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfabox_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfabox_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfabox params.input \
                                   -Vtk.OutputName test_md_embedded1d3d_1p_richards_tpfabox")

dune_add_test(NAME test_md_embedded1d3d_1p_richards_tpfabox_anadiff
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilBox
              COMPILE_DEFINITIONS DIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfabox_anadiff_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfabox_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfabox_anadiff_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_tpfabox_anadiff params.input \
                                   -Vtk.OutputName test_md_embedded1d3d_1p_richards_tpfabox_anadiff")

dune_add_test(NAME test_md_embedded1d3d_1p_richards_boxbox
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilBox
              COMPILE_DEFINITIONS ROOTTYPETAG=RootBox
              COMPILE_DEFINITIONS DIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_boxbox_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_boxbox_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_boxbox_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_boxbox params.input \
                                   -Vtk.OutputName test_md_embedded1d3d_1p_richards_boxbox")

dune_add_test(NAME test_md_embedded1d3d_1p_richards_boxbox_anadiff
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilBox
              COMPILE_DEFINITIONS ROOTTYPETAG=RootBox
              COMPILE_DEFINITIONS DIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_boxbox_anadiff_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_boxbox_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_boxbox_anadiff_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded1d3d_1p_richards_boxbox_anadiff params.input \
                                   -Vtk.OutputName test_md_embedded1d3d_1p_richards_boxbox_anadiff")
dune_symlink_to_source_files(FILES "params.input")
