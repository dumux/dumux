add_dumux_test(test_el2p test_el2p test_el2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/el2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/el2p-00003.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_el2p
  ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_el2p.dgf 100 200 100)

if(MPI_FOUND)
  add_dumux_test(test_el2p_parallel test_el2p_parallel test_el2p.cc
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/el2p-parallel-reference.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/s0004-p0001-el2p-00003.vtu
    "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_el2p_parallel"
    ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_el2p.dgf 100 200 100)
endif(MPI_FOUND)
