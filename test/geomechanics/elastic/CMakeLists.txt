add_input_file_links()

add_dumux_test(test_elastic test_elastic test_elastic.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/elasticmatrix-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/elasticmatrix-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_elastic"
                 --zeroThreshold {"uy":1e-15,"uz":1e-15,"stress X":1e-8,"stress Y":1e-8,"stress Z":1e-8})

#install sources
install(FILES
elasticmatrixproblem.hh
elasticspatialparams.hh
test_elastic.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/geomechanics/elastic)
