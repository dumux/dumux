add_dumux_test(test_general_box test_generalproblem2p test_generalproblem2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/generallens_box-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/generallens_box-00003.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_generalproblem2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p.input
  -ModelType Box
  -TimeManager.TEnd 1e2
  -TimeManager.DtInitial 2e1)

add_dumux_test(test_general_dec test_generalproblem2p test_generalproblem2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/generallens_decoupled-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/generallens_decoupled-00003.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_generalproblem2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p.input
  -ModelType Decoupled
  -TimeManager.TEnd 1e2
  -TimeManager.DtInitial 2e1)

