# the immiscible two-phase simulation program
#dune_add_test(NAME exercise1_2p_solution
#              SOURCES exercise1_2p_solution.cc
#              CMD_ARGS exercise1.input)

# the compositional two-phase two-component simulation program
dune_add_test(NAME exercise1_2p2c_solution
              SOURCES exercise1_2p2c.cc
              CMD_ARGS exercise1.input)

# the two-phase non-isothermal simulation program
dune_add_test(NAME exercise1_2pni_solution
              SOURCES exercise1_2pni_solution.cc
              CMD_ARGS exercise1.input)


# add tutorial to the common target
add_dependencies(build_tutorials exercise1_2p2c_solution exercise1_2pni_solution)

# add a symlink for the input file
dune_symlink_to_source_files(FILES "exercise1.input")
