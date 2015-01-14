if(MPI_CXX_FOUND)
  add_dumux_test(test_boxrichards test_boxrichards test_boxrichards.cc
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference-parallel.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-richardslensbox-00008.vtu
    "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards")
else()
# isothermal tests
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

# non-isothermal tests
add_dumux_test(test_boxrichardsniconvection test_boxrichardsniconvection test_boxrichardsniconvection.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectionbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection-00011.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection)

add_dumux_test(test_ccrichardsniconvection test_ccrichardsniconvection test_ccrichardsniconvection.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectioncc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconvection-00011.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconvection)

add_dumux_test(test_boxrichardsniconduction test_boxrichardsniconduction test_boxrichardsniconduction.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/richardsniconductionbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction)

add_dumux_test(test_ccrichardsniconduction test_ccrichardsniconduction test_ccrichardsniconduction.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/richardsniconductioncc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction)