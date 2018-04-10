add_input_file_links()

add_dumux_test(test_el2p test_el2p test_el2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/el2p-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/el2p-00047.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_el2p")

if(MPI_FOUND)
  add_dumux_test(test_el2p_parallel test_el2p test_el2p.cc
                 python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                   --script fuzzy
                   --files ${CMAKE_SOURCE_DIR}/test/references/el2p-parallel-reference.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/s0004-p0001-el2p-00047.vtu
                   --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_el2p")
endif()


# compile triaxial_el2p
add_dumux_test(triaxial_el2p triaxial_el2p triaxial_el2p.cc
  echo
  )

#install sources
install(FILES
el2pco2tables.hh
el2pproblem.hh
el2pspatialparams.hh
test_el2p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/geomechanics/el2p)
