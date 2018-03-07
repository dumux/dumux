dune_add_test(SOURCES test_fracture_2p_2p.cc)
dune_add_test(SOURCES noflowdomain_2p_ld.cc)
dune_symlink_to_source_files(FILES "grids" "test_fracture_2p_2p.input" "noflowdomain_2p_ld.input")
set(CMAKE_BUILD_TYPE Release)
