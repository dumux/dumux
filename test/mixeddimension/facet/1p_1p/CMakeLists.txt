dune_add_test(SOURCES test_fracture_1p_1p.cc)
dune_add_test(SOURCES test_analytical.cc)
dune_symlink_to_source_files(FILES "grids" "test_fracture_1p_1p.input" "test_analytical.input" "convergencetest.py")
set(CMAKE_BUILD_TYPE Release)
