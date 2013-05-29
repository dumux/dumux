# build target for the simple twophase lens problem in a general
# definition for both box and decoupled model
add_executable("test_generalproblem2p" EXCLUDE_FROM_ALL test_generalproblem2p.cc)

add_test(test_general_box
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/generallens_box-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/generallens_box-00003.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_generalproblem2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p.input
  -ModelType Box
  -TimeManager.TEnd 1e2
  -TimeManager.DtInitial 2e1)

add_test(test_general_dec
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/generallens_decoupled-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/generallens_decoupled-00003.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_generalproblem2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p.input
  -ModelType Decoupled
  -TimeManager.TEnd 1e2
  -TimeManager.DtInitial 2e1)

set(NORMALTESTS
  test_generalproblem2p)
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
