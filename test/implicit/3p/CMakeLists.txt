# isothermal tests
add_dumux_test(test_box3p test_box3p test_box3p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltration3pbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pbox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p)

add_dumux_test(test_cc3p test_cc3p test_cc3p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltration3pcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pcc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p)

# non-isothermal tests
add_dumux_test(test_box3pniconvection test_box3pniconvection test_box3pniconvection.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/3pniconvectionbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconvection-00011.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconvection)

add_dumux_test(test_cc3pniconvection test_cc3pniconvection test_cc3pniconvection.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/3pniconvectioncc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconvection-00011.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconvection)

add_dumux_test(test_box3pniconduction test_box3pniconduction test_box3pniconduction.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/3pniconductionbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconduction-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconduction)

add_dumux_test(test_cc3pniconduction test_cc3pniconduction test_cc3pniconduction.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/3pniconductioncc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconduction-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconduction)
