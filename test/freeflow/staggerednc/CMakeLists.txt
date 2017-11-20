add_input_file_links()

dune_add_test(NAME test_densitydrivenflow
              SOURCES test_densitydrivenflow.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-densitydriven-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_densitydrivenflow-00025.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_densitydrivenflow")

add_executable(test_stokes2c test_channel.cc)

dune_add_test(NAME test_purediffusion
              TARGET test_stokes2c
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-purediffusion.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_purediffusion-00013.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c test_purediffusion.input"
                             --zeroThreshold {"velocity_H2O \(m/s\)":1e-22})

dune_add_test(NAME test_advection
              TARGET test_stokes2c
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-advection.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_advection-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c test_advection.input")

add_executable(test_stokes2cni test_channel.cc)
target_compile_definitions(test_stokes2cni PUBLIC "NONISOTHERMAL=1")

dune_add_test(NAME test_stokes2cni_advection
              TARGET test_stokes2cni
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-advection.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni_advection-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}//test_stokes2cni test_stokes2cni_advection.input")

dune_add_test(NAME test_stokes2cni_diffusion
              TARGET test_stokes2cni
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-diffusion.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni_diffusion-00014.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}//test_stokes2cni test_stokes2cni_diffusion.input")
