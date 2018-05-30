dune_add_test(NAME test_1d3d_analytical2
              SOURCES test_1p_1p.cc)
dune_symlink_to_source_files(FILES "grids" "test_1p_1p.input")
set(CMAKE_BUILD_TYPE Release)
