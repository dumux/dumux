# build the tests for the bounding box tree
dumux_add_test(NAME test_bboxtree_dim1
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=1
              LABELS unit)

dumux_add_test(NAME test_bboxtree_dim2
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=2
              LABELS unit)

dumux_add_test(NAME test_bboxtree_dim3
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=3
              LABELS unit)

dumux_add_test(NAME test_geometry_fracture
              SOURCES test_geometry_fracture.cc
              CMAKE_GUARD dune-foamgrid_FOUND
              LABELS unit)

# symlink the input file in the build directory
dune_symlink_to_source_files(FILES "network1d.msh" "network2d.msh" "fracture.msh")
