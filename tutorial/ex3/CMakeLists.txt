# executables for exercise part a & b
# add_dumux_test(ex3_a ex3_a ex3_a.cc ex3_a)
# add_dumux_test(ex3_b ex3_b ex3_b.cc ex3_b)

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
