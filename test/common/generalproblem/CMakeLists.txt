add_dumux_test(test_general_box test_general_box test_generalproblem2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/generallens_box-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/generallens_box-00003.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_general_box
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p_reference.input
  -ModelType Box)

add_dumux_test(test_general_dec test_general_dec test_generalproblem2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/generallens_decoupled-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/generallens_decoupled-00003.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_general_dec
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p_reference.input
  -ModelType Decoupled)
