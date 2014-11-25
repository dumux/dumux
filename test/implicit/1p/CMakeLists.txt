add_dumux_test(test_box1p test_box1p test_box1p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box1p)

add_dumux_test(test_box1pwithamg test_box1pwithamg test_box1pwithamg.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/1ptestboxwithamg-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box1pwithamg)

add_dumux_test(test_cc1p test_cc1p test_cc1p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc1p)

add_dumux_test(test_cc1pwithamg test_cc1pwithamg test_cc1pwithamg.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/1ptestccwithamg-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc1pwithamg)
