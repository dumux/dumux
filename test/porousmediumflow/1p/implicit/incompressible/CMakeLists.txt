dune_symlink_to_source_files(FILES "test_1pincompressible.input")

# incompressible
add_dumux_test(test_1pincompressible test_1pincompressible test_1pincompressible.cc
               ${CMAKE_CURRENT_BINARY_DIR}/test_1pincompressible)

set(CMAKE_BUILD_TYPE Release)