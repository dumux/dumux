add_dumux_test(test_box2pni test_box2pni test_box2pni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injection2pnibox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2pni)

add_dumux_test(test_cc2pni test_cc2pni test_cc2pni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injection2pnicc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2pni)
