# the immiscible two-phase simulation program
add_dumux_test(exercise1_2p exercise1_2p exercise1_2p.cc ./exercise1_2p exercise1.input)

# the compositional two-phase simulation program
add_dumux_test(exercise1_2p2c exercise1_2p2c exercise1_2p2c.cc ./exercise1_2p2c exercise1.input)

# add tutorial to the common target
add_dependencies(build_tutorials exercise1_2p exercise1_2p2c)

# add a symlink for the input file
dune_symlink_to_source_files(FILES "exercise1.input")
