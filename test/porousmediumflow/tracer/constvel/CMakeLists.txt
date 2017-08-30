dune_symlink_to_source_files(FILES "test_tracer.input")

# tracer tests
dune_add_test(NAME test_tracer_tpfa
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa
              CMD_ARGS test_tracer.input -Problem.Name tracercc)

dune_add_test(NAME test_tracer_box
              SOURCES test_tracer.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox
              CMD_ARGS test_tracer.input -Problem.Name tracerbox)

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
tracertestproblem.hh
tracertestspatialparams.hh
test_tracer.cc
test_tracer.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/constvel)
