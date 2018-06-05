add_input_file_links()

dune_add_test(NAME test_ff_rans_bfs_kepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=BFSKEpsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_bfs_kepsilon-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_kepsilon-00070.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_kepsilon params.input
                        -Problem.Name test_ff_rans_bfs_kepsilon")


dune_add_test(NAME test_ff_rans_bfs_lowrekepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=BFSLowReKEpsilon
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_bfs_lowrekepsilon-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_lowrekepsilon-00053.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_lowrekepsilon params.input
                        -Problem.Name test_ff_rans_bfs_lowrekepsilon")

dune_add_test(NAME test_ff_rans_bfs_komega
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=BFSKOmega
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_bfs_komega-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_komega-00049.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_komega params.input
                        -Problem.Name test_ff_rans_bfs_komega")

dune_add_test(NAME test_ff_rans_bfs_oneeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=BFSOneEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_bfs_oneeq-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_oneeq-00044.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_oneeq params.input
                        -Problem.Name test_ff_rans_bfs_oneeq")

dune_add_test(NAME test_ff_rans_bfs_zeroeq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=BFSZeroEq
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_bfs_zeroeq-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_zeroeq-00035.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_bfs_zeroeq params.input
                        -Problem.Name test_ff_rans_bfs_zeroeq")
