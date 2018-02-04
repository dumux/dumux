dune_symlink_to_source_files(FILES "singlefracture.input" "single_equi.msh")
dune_add_test(NAME singlefracture_mpfa SOURCES test_singlefracture.cc)
set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
1ptestproblem.hh
1ptestspatialparams.hh
tracertestproblem.hh
tracertestspatialparams.hh
test_singlefracture.cc
singlefracture.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/1ptracer/singlefracture)
