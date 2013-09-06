# add build targets
add_executable("test_adaptive2p2c" EXCLUDE_FROM_ALL test_adaptive2p2c.cc)
add_dune_alugrid_flags(test_adaptive2p2c)
add_test(test_adaptive2p2c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_adaptive2p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_adaptive2p2c.input)

add_executable("test_dec2p2c" EXCLUDE_FROM_ALL test_dec2p2c.cc)
add_test(test_dec2p2c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_dec2p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c-00021.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_dec2p2c.input)

add_executable("test_multiphysics2p2c" EXCLUDE_FROM_ALL test_multiphysics2p2c.cc)
add_test(test_multiphysics2p2c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_multiphysics2p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c-00021.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c)

set(NORMALTESTS
  test_adaptive2p2c
  test_dec2p2c
  test_multiphysics2p2c)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY TEST_INCLUDE_FILE ${CMAKE_CURRENT_BINARY_DIR}/BuildTests.cmake)

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
