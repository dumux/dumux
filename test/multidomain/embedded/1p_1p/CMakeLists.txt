dune_add_test(SOURCES test_1p_1p.cc)
dune_symlink_to_source_files(FILES "grids" "test_1p_1p.input")
set(CMAKE_BUILD_TYPE Release)
