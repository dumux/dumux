dune_add_test(SOURCES test_1p2c_richards2c.cc)
dune_symlink_to_source_files(FILES "lupine.dgf" "test_1p2c_richards2c.input")
set(CMAKE_BUILD_TYPE Release)
