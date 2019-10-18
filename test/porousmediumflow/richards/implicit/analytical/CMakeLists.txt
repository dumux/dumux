add_input_file_links(FILES params.input)

# comparison to analytical solution - only with tpfa
dumux_add_test(NAME test_richards_analytical_tpfa
              LABELS porousmediumflow richards
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsAnalyticalCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richards_analytical_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richards_analytical_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richards_analytical_tpfa params.input -Problem.Name test_richards_analytical_tpfa")
