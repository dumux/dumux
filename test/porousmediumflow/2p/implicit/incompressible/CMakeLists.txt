dune_symlink_to_source_files(FILES "test_2pincompressible.input")

# incompressible
add_dumux_test(test_2pincompressible test_2pincompressible test_2pincompressible.cc
               ${CMAKE_CURRENT_BINARY_DIR}/test_2pincompressible)

set(CMAKE_BUILD_TYPE Release)