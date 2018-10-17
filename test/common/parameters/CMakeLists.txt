dune_add_test(SOURCES test_loggingparametertree.cc
              LABELS unit)
dune_symlink_to_source_files(FILES "params.input")

#install sources
install(FILES
test_loggingparametertree.cc
params.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/parameters)
