dune_add_test(SOURCES test_fracture_1p2c_1p2c.cc)
dune_symlink_to_source_files(FILES "grids" "test_fracture_1p2c_1p2c.input")
set(CMAKE_BUILD_TYPE Release)
