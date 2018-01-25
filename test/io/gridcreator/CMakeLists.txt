add_input_file_links()

dune_add_test(NAME test_gridcreator_gmsh
              SOURCES test_gridcreator_gmsh.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh"
                       --files ${CMAKE_SOURCE_DIR}/test/references/bifurcation-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/bifurcation-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation-00001.vtu)

dune_add_test(NAME test_gridcreator_cake
              SOURCES test_gridcreator_cake.cc
              CMAKE_GUARD "( dune-uggrid_FOUND OR dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_cake"
                       --files ${CMAKE_SOURCE_DIR}/test/references/cake-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/cake-00000.vtu)

dune_add_test(NAME test_gridcreator_subgrid
              SOURCES test_gridcreator_subgrid.cc
              CMAKE_GUARD dune-subgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_subgrid"
                       --files ${CMAKE_SOURCE_DIR}/test/references/subgrid-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/subgrid_three.vtu)
