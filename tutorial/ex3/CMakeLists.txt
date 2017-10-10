# executables for exercise part a & b
# dune_add_test(NAME ex3_a SOURCES ex3_a.cc)
# dune_add_test(NAME ex3_b SOURCES ex3_b.cc)

# add a symlink for the input file
dune_symlink_to_source_files(FILES "ex3_a.input" "ex3_b.input")

#install sources
install(FILES
ex3_a.cc
ex3_a_problem.hh
ex3_b.cc
ex3_b_problem.hh
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/tutorial/ex3)
