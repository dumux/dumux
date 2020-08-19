dune_symlink_to_source_files(FILES "params.input")

# executable for tpfa tests
add_executable(test_1p_convergence_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1p_convergence_tpfa PUBLIC "TYPETAG=OnePIncompressibleTpfa")

# executable for box tests
add_executable(test_1p_convergence_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1p_convergence_box PUBLIC "TYPETAG=OnePIncompressibleBox")

# using tpfa and conforming refinement
dumux_add_test(NAME test_1p_convergence_tpfa_conforming
              TARGET test_1p_convergence_tpfa
              CMAKE_GUARD "( dune-functions_FOUND )"
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_convergence_tpfa
              CMD_ARGS params.input -Problem.Name test_1p_convergence_tpfa_conforming)

# using tpfa and nonconforming refinement
dumux_add_test(NAME test_1p_convergence_tpfa_nonconforming
              TARGET test_1p_convergence_tpfa
              CMAKE_GUARD "( dune-functions_FOUND )"
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_convergence_tpfa
              CMD_ARGS params.input -Problem.Name test_1p_convergence_tpfa_conforming)

# using box and conforming refinement
dumux_add_test(NAME test_1p_convergence_box_conforming
              TARGET test_1p_convergence_box
              CMAKE_GUARD "( dune-functions_FOUND )"
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_convergence_box
              CMD_ARGS params.input -Problem.Name test_1p_convergence_box_conforming)

# using box and nonconforming refinement
dumux_add_test(NAME test_1p_convergence_box_nonconforming
              TARGET test_1p_convergence_box
              CMAKE_GUARD "( dune-functions_FOUND )"
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_convergence_box
              CMD_ARGS params.input -Problem.Name test_1p_convergence_box_conforming)
