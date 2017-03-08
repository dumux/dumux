dune_add_test(SOURCES test_fracture_2p_2p.cc)
dune_symlink_to_source_files(FILES "grids" "test_fracture_2p_2p.input")
set(CMAKE_BUILD_TYPE Release)
