add_subdirectory(1p_1p)
add_subdirectory(1p_richards)
add_subdirectory(1p2c_richards2c)
add_subdirectory(root_soil_benchmark)

dune_symlink_to_source_files(FILES "grids")
