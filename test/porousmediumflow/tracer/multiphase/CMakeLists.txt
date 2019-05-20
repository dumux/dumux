dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME test_tracer_multiphase_tpfa
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_multiphase_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_multiphase_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_multiphase_tpfa params.input -Problem.Name test_tracer_multiphase_tpfa")

dune_add_test(NAME test_tracer_multiphase_mpfa
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestMpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_multiphase_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_multiphase_mpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_multiphase_mpfa params.input -Problem.Name test_tracer_multiphase_mpfa")
