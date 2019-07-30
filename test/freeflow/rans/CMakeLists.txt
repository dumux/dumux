add_input_file_links()
dune_symlink_to_source_files(FILES references scripts)

add_executable(test_ff_rans_lauferpipe_kepsilon EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_rans_lauferpipe_kepsilon PUBLIC "TYPETAG=PipeLauferKEpsilon")

dune_add_test(NAME test_ff_rans_lauferpipe_kepsilon_twolayers
              TARGET test_ff_rans_lauferpipe_kepsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              TIMEOUT 1500
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_kepsilon_twolayers-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_kepsilon_twolayers-00070.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_kepsilon params.input
                        -Problem.Name test_ff_rans_lauferpipe_kepsilon_twolayers")

dune_add_test(NAME test_ff_rans_lauferpipe_kepsilon_wallfunction
              TARGET test_ff_rans_lauferpipe_kepsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              TIMEOUT 1500
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_kepsilon_wallfunction-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_kepsilon_wallfunction-00034.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_kepsilon params_wallfunction.input
                        -Problem.Name test_ff_rans_lauferpipe_kepsilon_wallfunction")

dune_add_test(NAME test_ff_rans_lauferpipe_lowrekepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=PipeLauferLowReKEpsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_lowrekepsilon-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_lowrekepsilon-00053.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_lowrekepsilon params.input
                        -Problem.Name test_ff_rans_lauferpipe_lowrekepsilon")

dune_add_test(NAME test_ff_rans_lauferpipe_komega
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=PipeLauferKOmega
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_komega-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_komega-00047.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_komega params.input
                        -Problem.Name test_ff_rans_lauferpipe_komega")

dune_add_test(NAME test_ff_rans_lauferpipe_oneeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=PipeLauferOneEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_oneeq-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_oneeq-00044.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_oneeq params.input
                        -Problem.Name test_ff_rans_lauferpipe_oneeq")

dune_add_test(NAME test_ff_rans_lauferpipe_zeroeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=PipeLauferZeroEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_zeroeq-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_zeroeq-00035.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_zeroeq params.input
                        -Problem.Name test_ff_rans_lauferpipe_zeroeq")

#///////////////////////////////
#/// The nonisothermal tests ///
#///////////////////////////////

add_executable(test_ff_ransni_lauferpipe_kepsilon EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_ransni_lauferpipe_kepsilon PUBLIC "TYPETAG=PipeLauferNIKEpsilon" "NONISOTHERMAL=1")

dune_add_test(NAME test_ff_ransni_lauferpipe_kepsilon_twolayers
              TARGET test_ff_ransni_lauferpipe_kepsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              TIMEOUT 1500
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_ransni_lauferpipe_kepsilon_twolayers-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_kepsilon_twolayers-00034.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_kepsilon params_nonisothermal.input
                        -Problem.Name test_ff_ransni_lauferpipe_kepsilon_twolayers")

dune_add_test(NAME test_ff_ransni_lauferpipe_kepsilon_wallfunction
              TARGET test_ff_ransni_lauferpipe_kepsilon
              LABELS freeflow rans
              TIMEOUT 1500
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_ransni_lauferpipe_kepsilon_wallfunction-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_kepsilon_wallfunction-00027.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_kepsilon params_nonisothermal_wallfunction.input
                        -Problem.Name test_ff_ransni_lauferpipe_kepsilon_wallfunction")

dune_add_test(NAME test_ff_ransni_lauferpipe_lowrekepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=PipeLauferNILowReKEpsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/./test_ff_ransni_lauferpipe_lowrekepsilon-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/./test_ff_ransni_lauferpipe_lowrekepsilon-00034.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_lowrekepsilon params_nonisothermal.input
                        -Problem.Name ./test_ff_ransni_lauferpipe_lowrekepsilon")
target_compile_definitions(test_ff_ransni_lauferpipe_lowrekepsilon PUBLIC "NONISOTHERMAL=1")

dune_add_test(NAME test_ff_ransni_lauferpipe_komega
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=PipeLauferNIKOmega
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_ransni_lauferpipe_komega-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_komega-00040.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_komega params_nonisothermal.input
                        -Problem.Name test_ff_ransni_lauferpipe_komega")
target_compile_definitions(test_ff_ransni_lauferpipe_komega PUBLIC "NONISOTHERMAL=1")

dune_add_test(NAME test_ff_ransni_lauferpipe_oneeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=PipeLauferNIOneEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_ransni_lauferpipe_oneeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_oneeq-00037.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_oneeq params_nonisothermal.input
                        -Problem.Name test_ff_ransni_lauferpipe_oneeq")
target_compile_definitions(test_ff_ransni_lauferpipe_oneeq PUBLIC "NONISOTHERMAL=1")

dune_add_test(NAME test_ff_ransni_lauferpipe_zeroeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=PipeLauferNIZeroEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_ransni_lauferpipe_zeroeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_zeroeq-00033.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_lauferpipe_zeroeq params_nonisothermal.input
                        -Problem.Name test_ff_ransni_lauferpipe_zeroeq")
target_compile_definitions(test_ff_ransni_lauferpipe_zeroeq PUBLIC "NONISOTHERMAL=1")
