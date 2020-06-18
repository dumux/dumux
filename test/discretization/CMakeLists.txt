add_subdirectory(cellcentered)
add_subdirectory(facecentered)
add_subdirectory(staggered)
add_subdirectory(box)
add_subdirectory(projection)
add_subdirectory(rotationsymmetry)

dumux_add_test(NAME test_disc_fvgridvariables
               LABELS unit
               SOURCES test_fvgridvariables.cc
               COMMAND ./test_disc_fvgridvariables
               CMD_ARGS -Problem.Name gridvarstest)
