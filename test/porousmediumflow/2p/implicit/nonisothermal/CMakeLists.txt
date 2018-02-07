dune_symlink_to_source_files(FILES "test_2pni.input")

dune_add_test(SOURCES test_2pni_fv.cc
              NAME test_2pni_box_simplex
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
              COMPILE_DEFINITIONS TYPETAG=InjectionBox2PNITypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-simplex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box_simplex-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box_simplex test_2pni.input -Problem.Name test_2pni_box_simplex -Grid.CellType Simplex")

dune_add_test(SOURCES test_2pni_fv.cc
              NAME test_2pni_box_cube
              COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
              COMPILE_DEFINITIONS TYPETAG=InjectionBox2PNITypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-cube-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box_cube-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box_cube test_2pni.input -Problem.Name test_2pni_box_cube")

dune_add_test(SOURCES test_2pni_fv.cc
              NAME test_2pni_tpfa_simplex
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
              COMPILE_DEFINITIONS TYPETAG=InjectionCC2PNITypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-simplex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa_simplex-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa_simplex test_2pni.input -Problem.Name test_2pni_tpfa_simplex -Grid.CellType Simplex")

dune_add_test(SOURCES test_2pni_fv.cc
              NAME test_2pni_tpfa_cube
              COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
              COMPILE_DEFINITIONS TYPETAG=InjectionCC2PNITypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-cube-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa_cube-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa_cube test_2pni.input -Problem.Name test_2pni_tpfa_cube")


set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem.hh
spatialparams.hh
test_2pni_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/nonisothermal)
