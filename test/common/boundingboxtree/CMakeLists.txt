# build the tests for the bounding box tree
dune_add_test(NAME test_bboxtree_dim1
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=1)

dune_add_test(NAME test_bboxtree_dim2
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=2)

dune_add_test(NAME test_bboxtree_dim3
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=3)

# symlink the input file in the build directory
dune_symlink_to_source_files(FILES "network1d.msh" "network2d.msh")

#install sources
install(FILES
test_bboxtree.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/boundingboxtree)

set(CMAKE_BUILD_TYPE Release)
