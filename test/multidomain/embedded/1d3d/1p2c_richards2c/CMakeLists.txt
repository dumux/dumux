dune_add_test(NAME test_embedded_1d3d_1p2c_richards2c
              SOURCES test_1p2c_richards2c.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UG )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p2c_richards2c_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_soilroot_1d-00003.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p2c_richards2c_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_soilroot_3d-00003.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_1d3d_1p2c_richards2c test_1p2c_richards2c.input")
dune_symlink_to_source_files(FILES "test_1p2c_richards2c.input")
set(CMAKE_BUILD_TYPE Release)
