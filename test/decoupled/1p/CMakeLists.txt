# add build targets
add_executable("test_dec1p" test_1p.cc)
add_test(test_dec1p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_1p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_1p-00001.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_1p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_1p.input)

add_executable("test_diffusion" test_diffusion.cc)
add_test(test_diffusion
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/diffusion-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/mimeticdiffusion-00001.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion
  3)

set(NORMALTESTS test_dec1p test_diffusion)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
