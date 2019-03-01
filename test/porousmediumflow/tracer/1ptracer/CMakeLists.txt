dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME test_1ptracer
              LABELS porousmediumflow tracer
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1ptracer_transport-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1ptracer-00010.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_1ptracer_pressure-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1p.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1ptracer params.input -Problem.Name test_1ptracer")
