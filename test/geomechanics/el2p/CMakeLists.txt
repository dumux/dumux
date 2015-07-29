add_input_file_links()

add_dumux_test(test_el2p test_el2p test_el2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/el2p-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/el2p-00003.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_el2p ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_el2p.dgf 100 200 100")

if(MPI_FOUND)
  add_dumux_test(test_el2p_parallel test_el2p_parallel test_el2p.cc
                 python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                   --script fuzzy
                   --files ${CMAKE_SOURCE_DIR}/test/references/el2p-parallel-reference.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/s0004-p0001-el2p-00003.vtu
                   --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_el2p_parallel ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_el2p.dgf 100 200 100")
endif()
