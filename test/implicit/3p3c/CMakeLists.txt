# isothermal tests
add_dumux_test(test_box3p3c test_box3p3c test_box3p3c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltrationbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/box3p3c-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p3c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p3c_reference.input)

add_dumux_test(test_cc3p3c test_cc3p3c test_cc3p3c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltrationcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/cc3p3c-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p3c_reference.input)

# non-isothermal tests
add_dumux_test(test_box3p3cnicolumnxylol test_box3p3cnicolumnxylol test_box3p3cnicolumnxylol.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/column3p3cnibox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/columnxylolbox-00058.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p3cnicolumnxylol)

add_dumux_test(test_cc3p3cnicolumnxylol test_cc3p3cnicolumnxylol test_cc3p3cnicolumnxylol.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/column3p3cnicc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/columnxylolcc-00054.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3cnicolumnxylol)

add_dumux_test(test_box3p3cnikuevette test_box3p3cnikuevette test_box3p3cnikuevette.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/kuevette3p3cnibox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/kuevettebox-00013.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p3cnikuevette)

# with Dune 2.3 we have 14 time steps, with Dune 2.4 we have 12
# remove this code if support for Dune 2.3 is dropped
set(REF_FILE_NUMBER "12")
if(("${DUNE_COMMON_VERSION_MAJOR}" STREQUAL "2") AND ("${DUNE_COMMON_VERSION_MINOR}" STREQUAL "3"))
  set(REF_FILE_NUMBER "14")
endif()
add_dumux_test(test_cc3p3cnikuevette test_cc3p3cnikuevette test_cc3p3cnikuevette.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/kuevette3p3cnicc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/kuevettecc-000${REF_FILE_NUMBER}.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3cnikuevette)
