dune_add_test(SOURCES test_analytical.cc)
dune_add_test(SOURCES heterogeneous.cc)
dune_symlink_to_source_files(FILES "grids" "test_analytical.input" "heterogeneous.input" "convergencetest.py")
set(CMAKE_BUILD_TYPE Release)
