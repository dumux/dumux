add_dumux_test(test_stokes2c test_stokes2c test_stokes2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2c-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c)

add_dumux_test(test_stokes2c_diffusion test_stokes2c test_stokes2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes2c_diffusion-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2c_diffusion-00020.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2c_diffusion.input)
