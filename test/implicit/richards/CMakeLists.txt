if(MPI_CXX_FOUND)
  add_dumux_test(test_boxrichards test_boxrichards test_boxrichards.cc
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference-parallel.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-richardslensbox-00008.vtu
    "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards")
else()
  add_dumux_test(test_boxrichards test_boxrichards test_boxrichards.cc
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/richardslensbox-00008.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards)
endif(MPI_CXX_FOUND)

add_dumux_test(test_ccrichards test_ccrichards test_ccrichards.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/richardslenscc-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichards)
