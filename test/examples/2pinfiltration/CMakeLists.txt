dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "initialsolutioncc.txt")

# using tpfa and point source
dune_add_test(NAME test_2p_pointsource_adaptive
              LABELS porousmediumflow 2p
              SOURCES main.cc
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_pointsource_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_pointsource_adaptive_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_pointsource_adaptive_tpfa params.input -Problem.Name test_2p_pointsource_adaptive_tpfa")
