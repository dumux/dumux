add_input_file_links()
dune_symlink_to_source_files(FILES laufer_re50000.csv laufer_re50000_u+y+.csv
                                   pdelab-kepsilon_twolayers.csv pdelab-kepsilon_wallfunction.csv
                                   pdelab-komega.csv pdelab-lowrekepsilon.csv pdelab-zeroeq.csv pdelab-oneeq.csv)

add_executable(test_pipe_laufer_kepsilon EXCLUDE_FROM_ALL test_pipe_laufer.cc)
target_compile_definitions(test_pipe_laufer_kepsilon PUBLIC "KEPSILON=1")

dune_add_test(NAME test_pipe_laufer_kepsilon_twolayers
              TARGET test_pipe_laufer_kepsilon
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_kepsilon.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_reference-00052.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer_kepsilon test_pipe_laufer_reference.input")

dune_add_test(NAME test_pipe_laufer_kepsilon_wallfunction
              TARGET test_pipe_laufer_kepsilon
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_kepsilon_wallfunction.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_reference-00029.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer_kepsilon test_pipe_laufer_reference_wallfunction.input")

dune_add_test(NAME test_pipe_laufer_lowrekepsilon
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_lowrekepsilon.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_reference-00066.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer_lowrekepsilon test_pipe_laufer_reference.input")
target_compile_definitions(test_pipe_laufer_lowrekepsilon PUBLIC "LOWREKEPSILON=1")

dune_add_test(NAME test_pipe_laufer_komega
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_komega.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_reference-00055.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer_komega test_pipe_laufer_reference.input")
target_compile_definitions(test_pipe_laufer_komega PUBLIC "KOMEGA=1")

dune_add_test(NAME test_pipe_laufer_oneeq
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_oneeq.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_reference-00046.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer_oneeq test_pipe_laufer_reference.input")
target_compile_definitions(test_pipe_laufer_oneeq PUBLIC "ONEEQ=1")

dune_add_test(NAME test_pipe_laufer_zeroeq
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_zeroeq.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_reference-00044.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer_zeroeq test_pipe_laufer_reference.input")

dune_add_test(NAME test_pipe_zeroeqni
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/pipe_zeroeqni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/pipe_zeroeqni-00034.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_zeroeqni test_pipe_zeroeqni.input")
target_compile_definitions(test_pipe_zeroeqni PUBLIC "NONISOTHERMAL=1")
