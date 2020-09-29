add_input_file_links()
dune_symlink_to_source_files(FILES solutions)

# isothermal tests
add_executable(test_ff_rans2c_kepsilon EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_rans2c_kepsilon PUBLIC "TYPETAG=FlatPlateNCKEpsilon")

dumux_add_test(NAME test_ff_rans2c_kepsilon_twolayer
              TARGET test_ff_rans2c_kepsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
                            COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2c_kepsilon_twolayer-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_kepsilon_twolayer-00032.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_kepsilon params.input
                             -Problem.Name test_ff_rans2c_kepsilon_twolayer")

dumux_add_test(NAME test_ff_rans2c_kepsilon_wallfunctions
              TARGET test_ff_rans2c_kepsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
                            COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2c_kepsilon_wallfunctions-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_kepsilon_wallfunctions-00028.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_kepsilon params_wallfunction.input
                             -Problem.Name test_ff_rans2c_kepsilon_wallfunctions")

dumux_add_test(NAME test_ff_rans2c_komega
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FlatPlateNCKOmega
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2c_komega-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_komega-00036.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_komega params.input
                             -Problem.Name test_ff_rans2c_komega")

dumux_add_test(NAME test_ff_rans2c_lowrekepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FlatPlateNCLowReKEpsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2c_lowrekepsilon-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_lowrekepsilon-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_lowrekepsilon params.input
                             -Problem.Name test_ff_rans2c_lowrekepsilon")

dumux_add_test(NAME test_ff_rans2c_oneeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FlatPlateNCOneEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2c_oneeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_oneeq-00032.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_oneeq params.input
                             -Problem.Name test_ff_rans2c_oneeq")


dumux_add_test(NAME test_ff_rans2c_zeroeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FlatPlateNCZeroEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2c_zeroeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_zeroeq-00030.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_zeroeq params.input
                             -Problem.Name test_ff_rans2c_zeroeq")

#///////////////////////////
#/// Nonisothermal tests ///
#///////////////////////////

add_executable(test_ff_rans2cni_kepsilon EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_rans2cni_kepsilon PUBLIC "TYPETAG=FlatPlateNCNIKEpsilon" "NONISOTHERMAL")

dumux_add_test(NAME test_ff_rans2cni_kepsilon_twolayer
              TARGET test_ff_rans2cni_kepsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_kepsilon_twolayer-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_kepsilon_twolayer-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_kepsilon params_nonisothermal.input
                             -Problem.Name test_ff_rans2cni_kepsilon_twolayer")

dumux_add_test(NAME test_ff_rans2cni_kepsilon_wallfunction
              TARGET test_ff_rans2cni_kepsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_kepsilon_wallfunction-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_kepsilon_wallfunction-00024.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_kepsilon params_nonisothermal_wallfunction.input
                             -Problem.Name test_ff_rans2cni_kepsilon_wallfunction")

dumux_add_test(NAME test_ff_rans2cni_komega
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FlatPlateNCNIKOmega
              COMPILE_DEFINITIONS NONISOTHERMAL
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_komega-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_komega-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_komega params_nonisothermal.input
                                       -Problem.Name test_ff_rans2cni_komega
                                       -Newton.TargetSteps 12")

dumux_add_test(NAME test_ff_rans2cni_komega_loadsolution
              TARGET test_ff_rans2cni_komega
              LABELS freeflow rans loadsolution
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_komega_loadsolution-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_komega_loadsolution-00008.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_komega params_nonisothermal.input
                                       -LoadSolution.RansCCFile ./solutions/komeganc_ccsolution.vtu
                                       -LoadSolution.RansFaceFile ./solutions/komeganc_facesolution.vtp
                                       -Problem.Name test_ff_rans2cni_komega_loadsolution
                                       -TimeLoop.DtInitial 5.0
                                       -TimeLoop.TEnd 200.0")

dumux_add_test(NAME test_ff_rans2cni_lowrekepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FlatPlateNCNILowReKEpsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_lowrekepsilon-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_lowrekepsilon-00030.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_lowrekepsilon params_nonisothermal.input
                             -Problem.Name test_ff_rans2cni_lowrekepsilon")
target_compile_definitions(test_ff_rans2cni_lowrekepsilon PUBLIC "NONISOTHERMAL")

dumux_add_test(NAME test_ff_rans2cni_oneeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FlatPlateNCNIOneEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_oneeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_oneeq-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_oneeq params_nonisothermal.input
                             -Problem.Name test_ff_rans2cni_oneeq")

target_compile_definitions(test_ff_rans2cni_oneeq PUBLIC "NONISOTHERMAL")

dumux_add_test(NAME test_ff_rans2cni_zeroeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FlatPlateNCNIZeroEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_zeroeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_zeroeq-00030.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_zeroeq params_nonisothermal.input
                             -Problem.Name test_ff_rans2cni_zeroeq")
target_compile_definitions(test_ff_rans2cni_zeroeq PUBLIC "NONISOTHERMAL")
