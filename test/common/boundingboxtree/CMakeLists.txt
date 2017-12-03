# build the tests for the bounding box tree
dune_add_test(SOURCES test_bboxtree.cc)

# symlink the input file in the build directory
dune_symlink_to_source_files(FILES "network1d.msh")

#install sources
install(FILES
test_bboxtree.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/boundingboxtree)
