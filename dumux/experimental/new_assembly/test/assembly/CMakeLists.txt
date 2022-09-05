add_subdirectory(poisson)

dumux_add_test(SOURCES test_matrix_pattern.cc LABELS unit new_assembly)
dumux_add_test(SOURCES test_residualfunction.cc LABELS unit new_assembly)
