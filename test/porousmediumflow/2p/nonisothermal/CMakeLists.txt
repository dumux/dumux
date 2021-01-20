dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_2pni_box_simplex
              LABELS porousmediumflow 2p  2pni
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
              COMPILE_DEFINITIONS TYPETAG=InjectionBox2PNITypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pni_box_simplex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box_simplex-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box_simplex params.input -Problem.Name test_2pni_box_simplex -Grid.CellType Simplex")

dumux_add_test(NAME test_2pni_box_cube
              LABELS porousmediumflow 2p  2pni
              SOURCES main.cc
              COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
              COMPILE_DEFINITIONS TYPETAG=InjectionBox2PNITypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pni_box_cube-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box_cube-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box_cube params.input -Problem.Name test_2pni_box_cube")

dumux_add_test(NAME test_2pni_tpfa_simplex
              LABELS porousmediumflow 2p  2pni
              SOURCES main.cc
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
              COMPILE_DEFINITIONS TYPETAG=InjectionCC2PNITypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pni_tpfa_simplex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa_simplex-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa_simplex params.input -Problem.Name test_2pni_tpfa_simplex -Grid.CellType Simplex")

dumux_add_test(NAME test_2pni_tpfa_cube
              LABELS porousmediumflow 2p  2pni
              SOURCES main.cc
              COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
              COMPILE_DEFINITIONS TYPETAG=InjectionCC2PNITypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pni_tpfa_cube-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa_cube-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa_cube params.input -Problem.Name test_2pni_tpfa_cube")
