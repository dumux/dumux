dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME example_1ptracer
               LABELS porousmediumflow tracer example
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1ptracer_transport-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/example_1ptracer-00010.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_1ptracer_pressure-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/example_1ptracer params.input -Problem.Name example_1ptracer")
