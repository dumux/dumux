dune_symlink_to_source_files(FILES "test_2pni.input")

# non-isothermal tests
if(HAVE_UG)
  dune_add_test(SOURCES test_2pni_box.cc
                NAME test_2pni_box
                COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                               --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-simplex-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box-00007.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box test_2pni.input -Problem.Name test_2pni_box -Grid.CellType Simplex")
else()
dune_add_test(SOURCES test_2pni_box.cc
              NAME test_2pni_box
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-cube-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box-00007.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_box test_2pni.input -Problem.Name test_2pni_box")
endif()

if(HAVE_UG)
dune_add_test(SOURCES test_2pni_cc.cc
              NAME test_2pni_tpfa
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/2pni_tpfa-00008.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa test_2pni.input -Problem.Name 2pni_tpfa -Grid.CellType Simplex")
else()
dune_add_test(SOURCES test_2pni_cc.cc
              NAME test_2pni_tpfa
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-cube-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/2pni_tpfa-00008.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pni_tpfa test_2pni.input -Problem.Name 2pni_tpfa")
endif()

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem.hh
spatialparams.hh
test_2pni_box.cc
test_2pni_cc.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/nonisothermal)
