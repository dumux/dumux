add_input_file_links()

add_dumux_test(test_stokes2c test_stokes2c test_stokes2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2c-00006.vtu
                  --zeroThreshold {"v":1e-9})

add_dumux_test(test_stokes2c_diffusion_mole test_stokes2c_mole test_stokes2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c_mole -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2c_diffusion.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c_diffusion-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2c_diffusion-00019.vtu
                  --zeroThreshold {"v":1e-16})

add_dumux_test(test_stokes2c_diffusion_mass test_stokes2c_mass test_stokes2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c_mass -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2c_diffusion.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c_diffusion-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2c_diffusion-00019.vtu
                  --zeroThreshold {"v":1e-16})
target_compile_definitions(test_stokes2c_mass PUBLIC "USE_MASSES=1")

#install sources
install(FILES
stokes2ctestproblem.hh
test_stokes2c.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/stokes2c)
