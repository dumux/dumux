dune_symlink_to_source_files(FILES "el1p.input")
dune_add_test(NAME test_el1p
              SOURCES test_el1p.cc
              COMMAND ./test_el1p
              CMD_ARGS el1p.input)
set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
1pproblem.hh
1pspatialparams.hh
poroelasticproblem.hh
poroelasticspatialparams.hh
test_el1p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/poromechanics/el1p)
