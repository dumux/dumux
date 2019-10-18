dune_symlink_to_source_files(FILES "params_advection.input" "params_diffusion.input" "params_advectionni.input" "params_diffusionni.input")

add_executable(test_ff_stokes2c_mass EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_stokes2c_mass PUBLIC "USE_MASS=1")

dumux_add_test(NAME test_ff_stokes2c_diffusion_mass
              TARGET test_ff_stokes2c_mass
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_diffusion-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_diffusion_mass-00015.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_mass params_diffusion.input
                             -Problem.Name test_ff_stokes2c_diffusion_mass"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-22})

add_executable(test_ff_stokes2c EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_ff_stokes2c_diffusion_mole
              TARGET test_ff_stokes2c
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_diffusion-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_diffusion_mole-00015.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c params_diffusion.input
                             -Problem.Name test_ff_stokes2c_diffusion_mole"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-22})

dumux_add_test(NAME test_ff_stokes2c_advection
              TARGET test_ff_stokes2c
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_advection-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_advection-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c params_advection.input
                             -Problem.Name test_ff_stokes2c_advection")

dumux_add_test(NAME test_ff_stokes2c_advection_nocaching
              SOURCES main.cc
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_advection-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_advection_nocaching-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_advection_nocaching params_advection.input
                             -Problem.Name test_ff_stokes2c_advection_nocaching")

add_executable(test_ff_stokes2cni EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_stokes2cni PUBLIC "NONISOTHERMAL=1")

dumux_add_test(NAME test_ff_stokes2cni_advection
              TARGET test_ff_stokes2cni
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2cni_advection-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2cni_advection-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2cni params_advectionni.input
                             -Problem.Name test_ff_stokes2cni_advection")

dumux_add_test(NAME test_ff_stokes2cni_diffusion
              TARGET test_ff_stokes2cni
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2cni_diffusion-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2cni_diffusion-00014.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2cni params_diffusionni.input
                             -Problem.Name test_ff_stokes2cni_diffusion"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-13})
