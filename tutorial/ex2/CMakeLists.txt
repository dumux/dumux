# the compositional two-phase simulation program
dune_add_test(NAME exercise2
              SOURCES exercise2.cc
              COMPILE_DEFINITIONS TYPETAG=Injection2p2pcCCProblem
              CMD_ARGS exercise2.input)

# add tutorial to the common target
add_dependencies(build_tutorials exercise2)

# add a symlink for the input file
dune_symlink_to_source_files(FILES "exercise2.input")
