dune_add_test(SOURCES test_fracture_1p2c_1p2c.cc)
dune_add_test(SOURCES test_fracture_1p2cni_1p2cni.cc)
dune_symlink_to_source_files(FILES "grids" "test_fracture_1p2c_1p2c.input" "test_fracture_1p2cni_1p2cni.input")
set(CMAKE_BUILD_TYPE Release)
