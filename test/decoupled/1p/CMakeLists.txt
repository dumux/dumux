add_dumux_test(test_dec1p test_dec1p test_1p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_1p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_1p-00001.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_dec1p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_1p.input)

add_dumux_test(test_diffusion test_diffusion test_diffusion.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/diffusion-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/mimeticdiffusion-00001.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion
  3)

add_dumux_test(test_diffusion3d test_diffusion3d test_diffusion3d.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_diffusion3d-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion3d_mimetic-00001.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion3d
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_diffusion3d_reference.input)
