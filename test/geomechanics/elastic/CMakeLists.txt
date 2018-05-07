dune_symlink_to_source_files(FILES "elastic.input")

# using box and numeric differentiation
dune_add_test(NAME test_elastic_box
              SOURCES test_elastic.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox NUMDIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/elasticbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/elastic-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_elastic_box elastic.input")

set(CMAKE_BUILD_TYPE Release)

install(FILES
problem.hh
spatialparams.hh
test_elastic.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/test/geomechanics/elastic)
