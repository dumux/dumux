# test: test_2cnistokes2p2cni_1ff
add_dumux_test(test_2cnistokes2p2cni_1ff test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cni-ff-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00004.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_reference.input
  )

# test: test_2cnistokes2p2cni_2pm
add_dumux_test(test_2cnistokes2p2cni_2pm test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cni-pm-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/2p2cni-00004.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_reference.input
  )

# test: test_2cnistokes2p2cni_boundarylayer_1ff
add_dumux_test(test_2cnistokes2p2cni_boundarylayer_1ff test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cniboundarylayer-ff-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00004.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_boundarylayer.input
  )

# test: test_2cnistokes2p2cni_boundarylayer_2pm
add_dumux_test(test_2cnistokes2p2cni_boundarylayer_2pm test_2cnistokes2p2cni test_2cnistokes2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2cnistokes2p2cniboundarylayer-pm-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/2p2cni-00004.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni_boundarylayer.input
  )
