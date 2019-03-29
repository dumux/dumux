add_input_file_links()

add_dumux_test(test_elrichards test_elrichards test_elrichards.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/elrichards-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/elrichards-00047.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_elrichards")

if(MPI_FOUND)
  add_dumux_test(test_elrichards_parallel test_elrichards test_elrichards.cc
                 python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                   --script fuzzy
                   --files ${CMAKE_SOURCE_DIR}/test/references/elrichards-parallel-reference.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/s0004-p0001-elrichards-00047.vtu
                   --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_elrichards")
endif()

#install sources
install(FILES
elrichardsco2tables.hh
elrichardsproblem.hh
elrichardsspatialparams.hh
test_elrichards.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/geomechanics/elrichards)
