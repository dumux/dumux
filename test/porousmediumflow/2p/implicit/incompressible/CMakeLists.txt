dune_symlink_to_source_files(FILES "test_2p.input")

# using tpfa
dune_add_test(NAME test_2p_incompressible_tpfa
              SOURCES test_2p_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleTpfa
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/lenscc-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/2p_tpfa-00008.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa test_2p.input -Problem.Name 2p_tpfa")

# using box
dune_add_test(NAME test_2p_incompressible_box
              SOURCES test_2p_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleBox
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/lensbox-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/2p_box-00007.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box test_2p.input -Problem.Name 2p_box")

# using mpfa
dune_add_test(NAME test_2p_incompressible_mpfa
              SOURCES test_2p_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleMpfa
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/lenscc-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/2p_mpfa-00008.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_mpfa test_2p.input -Problem.Name 2p_mpfa")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem.hh
spatialparams.hh
test_2p_box.cc
test_2p_cc.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/incompressible)
