dune_symlink_to_source_files(FILES "test_1pcompressible.input")

# incompressible
add_dumux_test(test_1pcompressible test_1pcompressible test_1pcompressible.cc
               ${CMAKE_CURRENT_BINARY_DIR}/test_1pcompressible)
