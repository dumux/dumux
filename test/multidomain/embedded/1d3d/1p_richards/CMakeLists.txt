dune_add_test(SOURCES test_1p_richards.cc)
dune_symlink_to_source_files(FILES "lupine.dgf" "test_1p_richards.input")
set(CMAKE_BUILD_TYPE Release)
