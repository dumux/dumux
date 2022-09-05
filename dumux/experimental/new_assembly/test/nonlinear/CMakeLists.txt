dumux_add_test(SOURCES test_scalar_nonlinear_function_matrix_free.cc LABELS unit new_assembly)

if (${HAVE_EIGEN3})
    dumux_add_test(SOURCES test_scalar_nonlinear_function.cc
                LINK_LIBRARIES Eigen3::Eigen
                LABELS unit new_assembly)
else ()
    dumux_add_test(SOURCES test_scalar_nonlinear_function.cc
                   LABELS unit new_assembly)
endif ()
