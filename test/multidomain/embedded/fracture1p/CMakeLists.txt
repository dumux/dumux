dune_add_test(NAME test_embedded_fracture1p
              SOURCES test_1p_1p.cc)
dune_symlink_to_source_files(FILES "fracture.msh" "test_1p_1p.input")
set(CMAKE_BUILD_TYPE Release)
