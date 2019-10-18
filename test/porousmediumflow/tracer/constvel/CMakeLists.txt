dune_symlink_to_source_files(FILES "params.input")

# explicit tracer tests (mass fractions)
dumux_add_test(NAME test_tracer_explicit_tpfa
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_explicit_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa params.input -Problem.Name test_tracer_explicit_tpfa")

dumux_add_test(NAME test_tracer_explicit_box
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_explicit_box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box params.input -Problem.Name test_tracer_explicit_box")

# explicit tracer tests (mole fractions, should yield same result)
dumux_add_test(NAME test_tracer_explicit_tpfa_mol
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false USEMOLES=true
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_explicit_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa_mol-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa_mol params.input -Problem.Name test_tracer_explicit_tpfa_mol")

dumux_add_test(NAME test_tracer_explicit_box_mol
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false USEMOLES=true
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_explicit_box-reference.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box_mol-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box_mol params.input -Problem.Name test_tracer_explicit_box_mol")

# implicit tracer tests
dumux_add_test(NAME test_tracer_implicit_tpfa
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=true USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_implicit_tpfa-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_tpfa params.input -Problem.Name test_tracer_implicit_tpfa")

dumux_add_test(NAME test_tracer_implicit_box
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=true USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_implicit_box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_box-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_box params.input -Problem.Name test_tracer_implicit_box")
