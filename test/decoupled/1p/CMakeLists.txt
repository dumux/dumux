# add build targets
add_executable("test_dec1p" EXCLUDE_FROM_ALL test_1p.cc)
add_test(test_dec1p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_1p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_1p-00001.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_dec1p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_1p.input)

add_executable("test_diffusion" EXCLUDE_FROM_ALL test_diffusion.cc)
add_test(test_diffusion
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/diffusion-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/mimeticdiffusion-00001.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion
  3)

add_executable("test_diffusion3d" EXCLUDE_FROM_ALL test_diffusion3d.cc)
add_dune_alugrid_flags(test_diffusion3d)
add_dune_ug_flags(test_diffusion3d)
add_test(test_diffusion3d
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_diffusion3d-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion3d_mimetic-00001.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion3d
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_diffusion3d.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_diffusion3d_ug1.dgf)

set(NORMALTESTS
  test_dec1p
  test_diffusion
  test_diffusion3d)
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
