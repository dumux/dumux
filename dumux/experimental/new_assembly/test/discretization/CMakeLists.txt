add_subdirectory(ccmpfa-o)
add_subdirectory(cctpfa)

dumux_add_test(SOURCES test_cctpfa_element_solution.cc LABELS unit new_assembly)
