dune_symlink_to_source_files(FILES "params.input" references)

dune_add_test(NAME test_shallowwater
              SOURCES main.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater)
