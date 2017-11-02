# a compositional two-phase simulation program
add_dumux_test(exercise2 exercise2 exercise2.cc exercise2)

# add tutorial to the common target
add_dependencies(build_tutorials exercise2)

# add a symlink for the input file
dune_symlink_to_source_files(FILES "exercise2.input")
