dune_symlink_to_source_files(FILES "test_1p.input" "test_1p_stationary.input")

# compressible instationary
dune_add_test(NAME test_1p_compressible_tpfa
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleTpfa
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1p_tpfa-00010.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_tpfa test_1p.input -Problem.Name 1p_tpfa")

dune_add_test(NAME test_1p_compressible_mpfa
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleMpfa
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1p_mpfa-00010.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_mpfa test_1p.input -Problem.Name 1p_mpfa")

dune_add_test(NAME test_1p_compressible_box
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleBox
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1p_box-00010.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_box test_1p.input -Problem.Name 1p_box")

# compressible stationary
dune_add_test(NAME test_1p_compressible_stationary_tpfa
              SOURCES test_1p_stationary.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleTpfa
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1p_stationary_tpfa-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_stationary_tpfa test_1p_stationary.input -Problem.Name 1p_stationary_tpfa")

dune_add_test(NAME test_1p_compressible_stationary_mpfa
              SOURCES test_1p_stationary.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleMpfa
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1p_stationary_mpfa-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_stationary_mpfa test_1p_stationary.input -Problem.Name 1p_stationary_mpfa")

dune_add_test(NAME test_1p_compressible_stationary_box
              SOURCES test_1p_stationary.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleBox
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1p_stationary_box-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_compressible_stationary_box test_1p_stationary.input -Problem.Name 1p_stationary_box")

set(CMAKE_BUILD_TYPE Release)
