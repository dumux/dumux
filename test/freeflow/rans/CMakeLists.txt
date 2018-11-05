add_input_file_links()
dune_symlink_to_source_files(FILES laufer_re50000.csv laufer_re50000_u+y+.csv)

add_executable(test_ff_rans_lauferpipe_kepsilon EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_rans_lauferpipe_kepsilon PUBLIC "KEPSILON=1")

dune_add_test(NAME test_ff_rans_lauferpipe_kepsilon_twolayers
              TARGET test_ff_rans_lauferpipe_kepsilon
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_kepsilon_twolayers-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_kepsilon_twolayers-00052.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_kepsilon params.input
                        -Problem.Name test_ff_rans_lauferpipe_kepsilon_twolayers")

dune_add_test(NAME test_ff_rans_lauferpipe_kepsilon_wallfunction
              TARGET test_ff_rans_lauferpipe_kepsilon
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_kepsilon_wallfunction-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_kepsilon_wallfunction-00029.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_kepsilon params_wallfunction.input
                        -Problem.Name test_ff_rans_lauferpipe_kepsilon_wallfunction")

dune_add_test(NAME test_ff_rans_lauferpipe_lowrekepsilon
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_lowrekepsilon-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_lowrekepsilon-00066.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_lowrekepsilon params.input
                        -Problem.Name test_ff_rans_lauferpipe_lowrekepsilon")
target_compile_definitions(test_ff_rans_lauferpipe_lowrekepsilon PUBLIC "LOWREKEPSILON=1")

dune_add_test(NAME test_ff_rans_lauferpipe_komega
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_komega-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_komega-00055.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_komega params.input
                        -Problem.Name test_ff_rans_lauferpipe_komega")
target_compile_definitions(test_ff_rans_lauferpipe_komega PUBLIC "KOMEGA=1")

dune_add_test(NAME test_ff_rans_lauferpipe_oneeq
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_oneeq-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_oneeq-00046.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_oneeq params.input
                        -Problem.Name test_ff_rans_lauferpipe_oneeq")
target_compile_definitions(test_ff_rans_lauferpipe_oneeq PUBLIC "ONEEQ=1")

dune_add_test(NAME test_ff_rans_lauferpipe_zeroeq
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_lauferpipe_zeroeq-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_zeroeq-00044.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_lauferpipe_zeroeq params.input
                        -Problem.Name test_ff_rans_lauferpipe_zeroeq")

dune_add_test(NAME test_ff_ransni_zeroeq
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_ransni_zeroeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_zeroeq-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_ransni_zeroeq params_nonisothermal.input
                        -Problem.Name test_ff_ransni_zeroeq")
target_compile_definitions(test_ff_ransni_zeroeq PUBLIC "NONISOTHERMAL=1")
