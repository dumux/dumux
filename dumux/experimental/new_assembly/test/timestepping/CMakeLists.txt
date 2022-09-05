dumux_add_test(SOURCES test_multistage_variables.cc LABELS unit new_assembly)
dumux_add_test(SOURCES test_scalar_nonlinear_timedependent_function_matrix_free.cc LABELS unit new_assembly)
dumux_add_test(SOURCES test_md_scalar_nonlinear_timedependent_function_matrix_free.cc LABELS unit new_assembly)

if (${HAVE_EIGEN3})
    dumux_add_test(
        SOURCES test_heat_equation_tpfa.cc
        LABELS unit new_assembly
        LINK_LIBRARIES Eigen3::Eigen
        COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
        CMD_ARGS  --script fuzzy
        --files ${CMAKE_SOURCE_DIR}/test/references/test_heat_equation_imp_euler-reference.vtu
        ${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_default_backend_imp_euler.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/test_heat_equation_rk4-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_default_backend_rk4.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/test_heat_equation_exp_euler-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_eigen_backend_exp_euler.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/test_heat_equation_crank_nicholson-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_flat_backend_crank_nicholson.vtu
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_tpfa"
    )
else()
    dumux_add_test(
        SOURCES test_heat_equation_tpfa.cc
        LABELS unit new_assembly
        COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
        CMD_ARGS  --script fuzzy
                  --files ${CMAKE_SOURCE_DIR}/test/references/test_heat_equation_imp_euler-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_default_backend_imp_euler.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/test_heat_equation_rk4-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_default_backend_rk4.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/test_heat_equation_crank_nicholson-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_flat_backend_crank_nicholson.vtu
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_heat_equation_tpfa"
    )
endif ()
