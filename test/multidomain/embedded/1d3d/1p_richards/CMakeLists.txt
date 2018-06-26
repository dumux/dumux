dune_add_test(NAME test_embedded_1p_richards_tpfatpfa
              SOURCES test_1p_richards.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1p_richards_tpfa_1d.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1p_richards_tpfatpfa_1d-00014.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1p_richards_tpfa_3d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1p_richards_tpfatpfa_3d-00014.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1p_richards_tpfatpfa test_1p_richards.input \
                                   -Problem.Name test_embedded_1p_richards_tpfatpfa")

dune_add_test(NAME test_embedded_1p_richards_tpfabox
              SOURCES test_1p_richards.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1p_richards_tpfa_1d.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1p_richards_tpfabox_1d-00017.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1p_richards_box_3d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1p_richards_tpfabox_3d-00017.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1p_richards_tpfabox test_1p_richards.input \
                                   -Problem.Name test_embedded_1p_richards_tpfabox")
dune_symlink_to_source_files(FILES "lupine.dgf" "test_1p_richards.input")
set(CMAKE_BUILD_TYPE Release)
