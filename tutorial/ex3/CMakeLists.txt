# executables for exercise part a & b
# add_dumux_test(ex3_a ex3_a ex3_a.cc ex3_a)
# add_dumux_test(ex3_b ex3_b ex3_b.cc ex3_b)

# add symlinks for the input files
add_input_file_links()

#install sources
install(FILES
ex3_a.cc
ex3_a_problem.hh
ex3_b.cc
ex3_b_problem.hh
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/tutorial/ex3)
