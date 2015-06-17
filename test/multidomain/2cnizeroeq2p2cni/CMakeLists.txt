file(COPY evaporationRates.gp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# test: test_2cnizeroeq2p2cni_1ff
add_dumux_test(test_2cnizeroeq2p2cni_1ff test_2cnizeroeq2p2cni test_2cnizeroeq2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2cnizeroeq2p2cni-ff-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2cni-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cnizeroeq2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnizeroeq2p2cni_reference.input
  )

# test: test_2cnizeroeq2p2cni_2pm
add_dumux_test(test_2cnizeroeq2p2cni_2pm test_2cnizeroeq2p2cni test_2cnizeroeq2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2cnizeroeq2p2cni-pm-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/2p2cni-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cnizeroeq2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnizeroeq2p2cni_reference.input
  )
