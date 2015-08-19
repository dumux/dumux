add_input_file_links()

add_dumux_test(test_stokes2c test_stokes2c test_stokes2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2c-00007.vtu
                  --zeroThreshold {"v":1e-9})

add_dumux_test(test_stokes2c_diffusion test_stokes2c test_stokes2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2c_diffusion.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c_diffusion-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2c_diffusion-00020.vtu)
