# isothermal tests
add_dumux_test(test_box2p test_box2p test_box2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/lensbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/lensbox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2p)
  
add_dumux_test(test_boxadaptive2p test_boxadaptive2p test_boxadaptive2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/lensboxadaptive-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/lensboxadaptive-00015.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxadaptive2p)

add_dumux_test(test_cc2p test_cc2p test_cc2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/lenscc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/lenscc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2p)
  
add_dumux_test(test_ccadaptive2p test_ccadaptive2p test_ccadaptive2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/lensccadaptive-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/lensccadaptive-00018.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccadaptive2p)  

# non-isothermal tests

if(HAVE_UG)
add_dumux_test(test_box2pni test_box2pni test_box2pni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-simplex-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injection2pnibox-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2pni)
else()
add_dumux_test(test_box2pni test_box2pni test_box2pni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-cube-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injection2pnibox-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2pni)
endif(HAVE_UG)

if(HAVE_UG)
add_dumux_test(test_cc2pni test_cc2pni test_cc2pni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-simplex-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injection2pnicc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2pni)
else()
add_dumux_test(test_cc2pni test_cc2pni test_cc2pni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-cube-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injection2pnicc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2pni)
endif(HAVE_UG)
