add_subdirectory(cellcentered)
add_subdirectory(staggered)
add_subdirectory(box)
add_subdirectory(projection)
add_subdirectory(rotationsymmetry)

dune_add_test(NAME test_fvgridvariables
              LABELS unit
              SOURCES test_fvgridvariables.cc
              COMMAND ./test_fvgridvariables
              CMD_ARGS -Problem.Name gridvarstest)
