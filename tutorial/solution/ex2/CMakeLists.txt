# the compositional two-phase simulation program
dune_add_test(NAME exercise2_solution
              SOURCES exercise2_solution.cc
              CMD_ARGS exercise2.input -Problem.OnlyPlotMaterialLaws false)

# add tutorial to the common target
add_dependencies(build_tutorials exercise2_solution)

# add a symlink for the input file
dune_symlink_to_source_files(FILES "exercise2.input")
