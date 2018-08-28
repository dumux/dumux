dune_symlink_to_source_files(FILES "test_tracer.input")

# explicit tracer tests (mass fractions)
dune_add_test(NAME test_tracer_explicit_tpfa
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-constvel-explicit-cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/tracer_explicit_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa test_tracer.input -Problem.Name tracer_explicit_tpfa")

dune_add_test(NAME test_tracer_explicit_box
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-constvel-explicit-box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/tracer_explicit_box-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box test_tracer.input -Problem.Name tracer_explicit_box")

# explicit tracer tests (mole fractions, should yield same result)
dune_add_test(NAME test_tracer_explicit_tpfa_mol
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false USEMOLES=true
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-constvel-explicit-cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/tracer_explicit_tpfa_mol-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa test_tracer.input -Problem.Name tracer_explicit_tpfa_mol")

dune_add_test(NAME test_tracer_explicit_box_mol
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false USEMOLES=true
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-constvel-explicit-box-reference.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/tracer_explicit_box_mol-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box test_tracer.input -Problem.Name tracer_explicit_box_mol")

# implicit tracer tests
dune_add_test(NAME test_tracer_implicit_tpfa
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=true USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-constvel-implicit-cc-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/tracer_implicit_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_tpfa test_tracer.input -Problem.Name tracer_implicit_tpfa")

dune_add_test(NAME test_tracer_implicit_box
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=true USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-constvel-implicit-box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/tracer_implicit_box-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_box test_tracer.input -Problem.Name tracer_implicit_box")

#install sources
install(FILES
tracertestproblem.hh
tracertestspatialparams.hh
test_tracer.cc
test_tracer.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/constvel)
