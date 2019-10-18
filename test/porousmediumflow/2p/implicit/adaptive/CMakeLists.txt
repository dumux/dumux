dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "initialsolutioncc.txt")
dune_symlink_to_source_files(FILES "initialsolutionbox.txt")

# using tpfa
dumux_add_test(NAME test_2p_adaptive_tpfa
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveTpfa
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_tpfa params.input -Problem.Name test_2p_adaptive_tpfa")

# using tpfa and point source
dumux_add_test(NAME test_2p_pointsource_adaptive_tpfa
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPAdaptivePointSource
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_pointsource_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_pointsource_adaptive_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_pointsource_adaptive_tpfa params.input -Problem.Name test_2p_pointsource_adaptive_tpfa")

# using mpfa
dumux_add_test(NAME test_2p_adaptive_mpfa
              LABELS porousmediumflow 2p
             SOURCES main.cc
             COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveMpfa
             CMAKE_GUARD dune-alugrid_FOUND
             COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
             CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_mpfa-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_mpfa-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_mpfa params.input -Problem.Name test_2p_adaptive_mpfa")

# using box
dumux_add_test(NAME test_2p_adaptive_box
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveBox
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box params.input -Problem.Name test_2p_adaptive_box")
