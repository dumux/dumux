add_input_file_links()

dune_add_test(NAME test_stokes2c_densitydrivenflow
              SOURCES test_densitydrivenflow.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-densitydriven-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_densitydrivenflow-00021.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c_densitydrivenflow"
                             --zeroThreshold {"X^Air_liquid":1e-12}
                             --zeroThreshold {"x^Air_liquid":1e-12}
                             --zeroThreshold {"velocity_liquid \(m/s\)":1e-10}
                             --zeroThreshold {"deltaRho":1e-7})

add_executable(test_stokes2c_mass EXCLUDE_FROM_ALL test_channel.cc)
target_compile_definitions(test_stokes2c_mass PUBLIC "USE_MASS=1")

add_executable(test_stokes2c EXCLUDE_FROM_ALL test_channel.cc)

dune_add_test(NAME test_stokes2c_purediffusion_mass
              TARGET test_stokes2c_mass
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-purediffusion.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_purediffusion-00014.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c_mass test_stokes2c_purediffusion.input"
                             --zeroThreshold {"velocity_liquid \(m/s\)":1e-22})

dune_add_test(NAME test_stokes2c_purediffusion_mole
              TARGET test_stokes2c
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-purediffusion.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_purediffusion-00013.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c test_stokes2c_purediffusion.input"
                             --zeroThreshold {"velocity_liquid \(m/s\)":1e-22})

dune_add_test(NAME test_stokes2c_advection
              TARGET test_stokes2c
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-advection.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_advection-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c test_stokes2c_advection.input")

dune_add_test(NAME test_stokes2c_advection_no_caching
              SOURCES test_channel.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2c-advection.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_advection-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c_advection_no_caching test_stokes2c_advection.input")

add_executable(test_stokes2cni EXCLUDE_FROM_ALL test_channel.cc)
target_compile_definitions(test_stokes2cni PUBLIC "NONISOTHERMAL=1")

dune_add_test(NAME test_stokes2cni_advection
              TARGET test_stokes2cni
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-advection.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni_advection-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}//test_stokes2cni test_stokes2cni_advection.input")

dune_add_test(NAME test_stokes2cni_diffusion
              TARGET test_stokes2cni
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-diffusion.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni_diffusion-00014.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}//test_stokes2cni test_stokes2cni_diffusion.input"
                             --zeroThreshold {"velocity_liquid \(m/s\)":1e-13})

dune_add_test(NAME test_msfreeflow
              SOURCES test_msfreeflow.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/maxwellstefanfreeflow-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_msfreeflow-00005.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_msfreeflow test_msfreeflow.input")
