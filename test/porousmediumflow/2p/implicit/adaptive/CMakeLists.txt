dune_symlink_to_source_files(FILES "test_2p_adaptive.input")

# using tpfa
dune_add_test(NAME test_2p_adaptive_tpfa
              SOURCES test_2p_adaptive_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveTpfa
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/lensccadaptive-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/2p_adaptive_tpfa-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_tpfa test_2p_adaptive.input -Problem.Name 2p_adaptive_tpfa")

# using mpfa
dune_add_test(NAME test_2p_adaptive_mpfa
             SOURCES test_2p_adaptive_fv.cc
             COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveMpfa
             CMAKE_GUARD dune-alugrid_FOUND
             COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
             CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/lensccadaptive-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/2p_adaptive_mpfa-00008.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_mpfa test_2p_adaptive.input -Problem.Name 2p_adaptive_mpfa")

# using box
dune_add_test(NAME test_2p_adaptive_box
              SOURCES test_2p_adaptive_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveBox
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/lensboxadaptive-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/2p_adaptive_box-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box test_2p_adaptive.input -Problem.Name 2p_adaptive_box")

set(CMAKE_BUILD_TYPE Release)
