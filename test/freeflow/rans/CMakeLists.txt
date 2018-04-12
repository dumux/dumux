add_input_file_links()
dune_symlink_to_source_files(FILES laufer_re50000.csv laufer_re50000_u+y+.csv pdelab-lowrekepsilon.csv pdelab-zeroeq.csv)

dune_add_test(NAME test_pipe_laufer_zeroeq
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_zeroeq.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_reference-00044.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer_zeroeq test_pipe_laufer_reference.input")

dune_add_test(NAME test_pipe_laufer_lowrekepsilon
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_lowrekepsilon.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_reference-00074.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer_lowrekepsilon test_pipe_laufer_reference.input")
target_compile_definitions(test_pipe_laufer_lowrekepsilon PUBLIC "LOWREKEPSILON=1")

dune_add_test(NAME test_pipe_zeroeqni
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/pipe_zeroeqni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/pipe_zeroeqni-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_zeroeqni test_pipe_zeroeqni.input")
target_compile_definitions(test_pipe_zeroeqni PUBLIC "NONISOTHERMAL=1")
