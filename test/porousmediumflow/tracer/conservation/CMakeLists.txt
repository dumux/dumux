add_input_file_links()

add_executable(test_tracer_conservation EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_tracer_conservation_dirichlet
               TARGET test_tracer_conservation
               LABELS porousmediumflow tracer
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_conservation
               CMD_ARGS -Problem.UseDirichletForTracer true)

dumux_add_test(NAME test_tracer_conservation_neumann
               TARGET test_tracer_conservation
               LABELS porousmediumflow tracer
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_conservation
               CMD_ARGS -Problem.UseDirichletForTracer false)
