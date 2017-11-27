dune_symlink_to_source_files(FILES "test_tracer.input")

# tracer tests
dune_add_test(NAME test_tracer_explicit_tpfa
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false
              CMD_ARGS test_tracer.input -Problem.Name tracer_explicit_tpfa)

dune_add_test(NAME test_tracer_explicit_box
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false
              CMD_ARGS test_tracer.input -Problem.Name tracer_explicit_box)

dune_add_test(NAME test_tracer_implicit_tpfa
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=true
              CMD_ARGS test_tracer.input -Problem.Name tracer_implicit_tpfa)

dune_add_test(NAME test_tracer_implicit_box
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=true
              CMD_ARGS test_tracer.input -Problem.Name tracer_implicit_box)

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
tracertestproblem.hh
tracertestspatialparams.hh
test_tracer.cc
test_tracer.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/constvel)
