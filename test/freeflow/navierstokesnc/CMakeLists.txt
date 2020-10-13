add_subdirectory(densitydrivenflow)
add_subdirectory(channel)
add_subdirectory(maxwellstefan)

add_executable(test_ff_stokes2c_mass_new EXCLUDE_FROM_ALL main_new.cc)
target_compile_definitions(test_ff_stokes2c_mass_new PUBLIC "USE_MASS=1")

dumux_add_test(NAME test_ff_stokes2c_diffusion_mass_new
              TARGET test_ff_stokes2c_mass_new
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_diffusion-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_diffusion_mass_new-00015.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_mass_new channel/params_diffusion.input
                             -Problem.Name test_ff_stokes2c_diffusion_mass_new"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-22})

add_executable(test_ff_stokes2c_mole_new EXCLUDE_FROM_ALL main_new.cc)
target_compile_definitions(test_ff_stokes2c_mole_new PUBLIC "USE_MASS=0")

dumux_add_test(NAME test_ff_stokes2c_diffusion_mole_new
              TARGET test_ff_stokes2c_mole_new
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_diffusion-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_diffusion_mole_new-00015.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_mole_new channel/params_diffusion.input
                             -Problem.Name test_ff_stokes2c_diffusion_mole_new"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-22})

dumux_add_test(NAME test_ff_stokes2c_advection_new
              TARGET test_ff_stokes2c_mole_new
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_advection-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_advection_new-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_mole_new channel/params_advection.input
                             -Problem.Name test_ff_stokes2c_advection_new")

add_executable(test_ff_stokes2cni_new EXCLUDE_FROM_ALL main_new.cc)
target_compile_definitions(test_ff_stokes2cni_new PUBLIC "NONISOTHERMAL=1")

dumux_add_test(NAME test_ff_stokes2cni_advection_new
              TARGET test_ff_stokes2cni_new
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2cni_advection-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2cni_advection_new-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2cni_new channel/params_advectionni.input
                             -Problem.Name test_ff_stokes2cni_advection_new")

dumux_add_test(NAME test_ff_stokes2cni_diffusion_new
              TARGET test_ff_stokes2cni_new
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2cni_diffusion-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2cni_diffusion_new-00014.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2cni_new channel/params_diffusionni.input
                             -Problem.Name test_ff_stokes2cni_diffusion_new"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-13})

dumux_add_test(NAME test_ff_stokes2c_densitydrivenflow_new
           LABELS freeflow navierstokes
           TARGET test_ff_stokes2c_mole_new
           CMAKE_GUARD HAVE_UMFPACK
           COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
           CMD_ARGS       --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_densitydrivenflow-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_densitydrivenflow_new-00021.vtu
                          --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_mole_new densitydrivenflow/params.input
                          -Problem.Name test_ff_stokes2c_densitydrivenflow_new
                          -Problem.IsChannelTest false"
                          --zeroThreshold {"X^Air_liq":1e-12}
                          --zeroThreshold {"x^Air_liq":1e-12}
                          --zeroThreshold {"velocity_liq \(m/s\)":1e-10}
                          --zeroThreshold {"deltaRho":1e-7})
