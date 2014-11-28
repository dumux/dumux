# isothermal tests
add_dumux_test(test_box2p2c test_box2p2c test_box2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injectionbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injectionbox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2p2c)

add_dumux_test(test_cc2p2c test_cc2p2c test_cc2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injectioncc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injectioncc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2p2c)

# non-isothermal tests
add_dumux_test(test_box2p2cni test_box2p2cni test_box2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/waterairbox-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2p2cni)

add_dumux_test(test_cc2p2cni test_cc2p2cni test_cc2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/wateraircc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/wateraircc-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2p2cni)
