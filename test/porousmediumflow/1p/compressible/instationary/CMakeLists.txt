dune_symlink_to_source_files(FILES "params.input")

# compressible instationary
dumux_add_test(NAME test_1p_compressible_instationary_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_instationary_tpfa-00010.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_instationary_tpfa params.input -Problem.Name test_1p_compressible_instationary_tpfa")

dumux_add_test(NAME test_1p_compressible_instationary_mpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleMpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_instationary_mpfa-00010.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_instationary_mpfa params.input -Problem.Name test_1p_compressible_instationary_mpfa")

dumux_add_test(NAME test_1p_compressible_instationary_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_instationary_box-00010.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_instationary_box params.input -Problem.Name test_1p_compressible_instationary_box")
