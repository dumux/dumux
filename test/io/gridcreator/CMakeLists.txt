add_dumux_test(test_gridcreator_gmsh test_gridcreator_gmsh test_gridcreator_gmsh.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh"
                 --files ${CMAKE_SOURCE_DIR}/test/references/bifurcation-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/bifurcation-00000.vtu
                         ${CMAKE_SOURCE_DIR}/test/references/bifurcation-reference-refined.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/bifurcation-00001.vtu)
                         
add_dumux_test(test_gridcreator_cake_powerlaw test_gridcreator_cake test_gridcreator_cake.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_cake"
                 --files ${CMAKE_SOURCE_DIR}/test/references/cake-powerlaw-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cake-00000.vtu)

add_dumux_test(test_gridcreator_pieceofcake_lin test_gridcreator_cake test_gridcreator_cake.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_cake -Grid.PolarCoorMax '1000.0 10.0 1010.0' -Grid.UsePowerLawSpacerRadius false"
                 --files ${CMAKE_SOURCE_DIR}/test/references/pieceofcake-linspace-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/pieceofcake-linspace-00000.vtu)