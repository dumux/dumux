add_input_file_links()
dune_symlink_to_source_files(FILES convergence.sh)
set(CMAKE_BUILD_TYPE Release)
add_executable(test_closedsystem EXCLUDE_FROM_ALL test_closedsystem.cc)

dune_add_test(NAME test_navierstokes_liddrivencavity_re1
              TARGET test_closedsystem
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/liddrivencavity-low-re-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_liddrivencavity_re1-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_liddrivencavity_re1.input")

dune_add_test(NAME test_navierstokes_liddrivencavity_re1000
              TARGET test_closedsystem
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/liddrivencavity-high-re-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_liddrivencavity_re1000-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_liddrivencavity_re1000.input")

dune_add_test(NAME test_navierstokes_hydrostaticpressure
              TARGET test_closedsystem
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokeshydrostaticpressure-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_hydrostaticpressure-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_hydrostaticpressure.input"
                             --zeroThreshold {"velocity_Constant \(m/s\)":1e-16})

dune_add_test(NAME test_channel_stokes
              SOURCES test_channel.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/channel-stokes.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokes-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokes")

add_executable(test_channel_stokesni EXCLUDE_FROM_ALL test_channel.cc)
target_compile_definitions(test_channel_stokesni PUBLIC "NONISOTHERMAL=1")

dune_add_test(NAME test_channel_stokesni_convection
              TARGET test_channel_stokesni
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokesni-convection-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokesni_convection-00006.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokesni test_channel_stokesni_convection.input")

dune_add_test(NAME test_channel_stokesni_conduction
              TARGET test_channel_stokesni
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokesni-conduction-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokesni_conduction-00004.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokesni test_channel_stokesni_conduction.input"
                             --zeroThreshold {"velocity_H2O \(m/s\)":1e-20})

dune_add_test(NAME test_channel_navierstokes
              SOURCES test_channel.cc
              COMPILE_DEFINITIONS ENABLE_NAVIERSTOKES=1
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/channel-navierstokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel_navierstokes-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_navierstokes -Vtk.WriteFaceData 1")

dune_add_test(NAME test_channel_navierstokes_restart
              TARGET test_channel_navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/channel-navierstokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel_navierstokes_restart-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_navierstokes -Vtk.WriteFaceData 1 -TimeLoop.DtInitial 1 -Restart.Time 1 -CellCenter.Restart.File test_channel_navierstokes-00001.vtu -Face.Restart.File test_channel_navierstokes-face-00001.vtp  -Problem.Name test_channel_navierstokes_restart")

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_channel_navierstokes_restart PROPERTIES DEPENDS test_channel_navierstokes)

dune_add_test(NAME test_stokes_donea
              SOURCES test_donea.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes-donea-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_donea-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes_donea")

dune_add_test(NAME test_navierstokes_1d
              SOURCES test_navierstokes_1d.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_navierstokes_1d.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_navierstokes_1d-00001.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_navierstokes_1d")

dune_add_test(NAME test_navierstokes_kovasznay
              SOURCES test_kovasznay.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_kovasznay-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_kovasznay-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_navierstokes_kovasznay")

dune_add_test(NAME test_navierstokes_angeli
              SOURCES test_angeli.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_angeli-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_angeli-00045.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_navierstokes_angeli")

dune_add_test(NAME test_stokes_channel_3d
              SOURCES test_stokes_channel_3d.cc
              COMPILE_DEFINITIONS DIM_3D=1
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes_channel_3d-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes_channel_3d-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes_channel_3d"
                             --zeroThreshold {"velocity_component \(m/s\)":1e-12})

dune_add_test(NAME test_stokes_channel_pseudo3d
              SOURCES test_stokes_channel_3d.cc
              COMPILE_DEFINITIONS DIM_3D=0
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/stokes_channel_pseudo3d-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes_channel_pseudo3d-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes_channel_pseudo3d"
                             --zeroThreshold {"velocity_component \(m/s\)":1e-12})
