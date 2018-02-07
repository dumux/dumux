# executables for exercise part a & b
#part a: 2pproblem
dune_add_test(NAME exercise3_a
              SOURCES exercise3.cc
              COMPILE_DEFINITIONS TYPETAG=ExerciseThreeBoxTwoPTypeTag
              COMPILE_ONLY # for testing purposes, ignore for the exercise
              CMD_ARGS exercise3_a.input)

#part b: 2p2cproblem
dune_add_test(NAME exercise3_b
              SOURCES exercise3.cc
              COMPILE_DEFINITIONS TYPETAG=ExerciseThreeBoxTwoPTwoCTypeTag
              COMPILE_ONLY # for testing purposes, ignore for the exercise
              CMD_ARGS exercise3_b.input)

# add tutorial to the common target
add_dependencies(build_tutorials exercise3_a exercise3_b)

# add symlinks for the input files
add_input_file_links()
