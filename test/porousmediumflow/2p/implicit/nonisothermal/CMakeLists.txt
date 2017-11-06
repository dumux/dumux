dune_symlink_to_source_files(FILES "test_2pni.input")

# non-isothermal tests
if(HAVE_UG)
  dune_add_test(SOURCES test_box2pni.cc
                NAME test_box2pni
                COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                               --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-simplex-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/injection2pnibox-00007.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pni test_2pni.input -Problem.Name injection2pnibox -Grid.CellType Simplex")
else()
dune_add_test(SOURCES test_box2pni.cc
              NAME test_box2pni
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-cube-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/injection2pnibox-00007.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pni test_2pni.input -Problem.Name injection2pnibox")
endif()

if(HAVE_UG)
dune_add_test(SOURCES test_cc2pni.cc
              NAME test_cc2pni
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/injection2pnicc-00008.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pni test_2pni.input -Problem.Name injection2pnicc -Grid.CellType Simplex")
else()
dune_add_test(SOURCES test_cc2pni.cc
              NAME test_cc2pni
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-cube-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/injection2pnicc-00008.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pni test_2pni.input -Problem.Name injection2pnicc")
endif()

set(CMAKE_BUILD_TYPE Release)
