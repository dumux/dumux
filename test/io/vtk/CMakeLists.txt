dune_symlink_to_source_files(FILES "test-in.vtu")

set(CMAKE_BUILD_TYPE Debug)

dune_add_test(SOURCES test_vtkreader.cc
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader"
                       --files ${CMAKE_CURRENT_BINARY_DIR}/test-in.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test-out.vtu)
