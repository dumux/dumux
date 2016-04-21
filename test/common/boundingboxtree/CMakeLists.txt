# build the tests for the bounding box tree
add_dumux_test(test_bboxtree test_bboxtree test_bboxtree.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_bboxtree)

# symlink the input file in the build directory
dune_symlink_to_source_files(FILES "network1d.msh")

#install sources
install(FILES
test_bboxtree.cc
test_bboxtree.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/boundingboxtree)