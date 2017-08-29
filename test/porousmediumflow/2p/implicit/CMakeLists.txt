add_subdirectory(pointsources)
add_subdirectory(incompressible)
add_input_file_links()

# isothermal tests
add_dumux_test(test_box2p test_box2p test_box2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/lensbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/lensbox-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2p")

add_dumux_test(test_cc2p test_cc2p test_cc2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/lenscc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/lenscc-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2p")

add_dumux_test(test_ccadaptive2p test_ccadaptive2p test_ccadaptive2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/lensccadaptive-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/lensccadaptive-00017.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccadaptive2p")

add_dumux_test(test_boxadaptive2p test_boxadaptive2p test_boxadaptive2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/lensboxadaptive-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/lensboxadaptive-00014.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxadaptive2p")

add_dumux_test(test_generalizeddirichlet test_generalizeddirichlet test_generalizeddirichlet.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/generalizeddirichlet-reference.vtp
                         ${CMAKE_CURRENT_BINARY_DIR}/generalizeddirichlet-00034.vtp
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_generalizeddirichlet")

if(HAVE_OPM_GRID)
add_dumux_test(test_cc2pcornerpoint test_cc2pcornerpoint test_cc2pcornerpoint.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/cc2pcornerpoint-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc2pcornerpoint-00005.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pcornerpoint")
endif()

# non-isothermal tests
if(HAVE_UG)
add_dumux_test(test_box2pni test_box2pni test_box2pni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-simplex-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injection2pnibox-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pni -Grid.CellType Simplex")
else()
add_dumux_test(test_box2pni test_box2pni test_box2pni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-cube-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injection2pnibox-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pni")
endif()

if(HAVE_UG)
add_dumux_test(test_cc2pni test_cc2pni test_cc2pni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-simplex-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injection2pnicc-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pni -Grid.CellType Simplex")
else()
add_dumux_test(test_cc2pni test_cc2pni test_cc2pni.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-cube-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/injection2pnicc-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pni")
endif()

add_dumux_test(test_fracture_box2p test_fracture_box2p test_fracture_box2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/fracturebox2p-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/fracturebox-00023.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_fracture_box2p")

add_dumux_test(test_fracture_cc2p test_fracture_cc2p test_fracture_cc2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/fracturecc2p-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/fracturecc-00029.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_fracture_cc2p")

add_dumux_test(test_fracture_ccmpfa2p test_fracture_ccmpfa2p test_fracture_ccmpfa2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/fractureccmpfa2p-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/fractureccmpfa-00031.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_fracture_ccmpfa2p")

#install sources
install(FILES
cc2pcornerpointproblem.hh
cc2pcornerpointspatialparams.hh
generalizeddirichletproblem.hh
generalizeddirichletspatialparams.hh
injectionproblem2pni.hh
lensproblem.hh
lensspatialparams.hh
test_box2p.cc
test_box2pni.cc
test_boxadaptive2p.cc
test_cc2p.cc
test_cc2pcornerpoint.cc
test_cc2pni.cc
test_ccadaptive2p.cc
test_generalizeddirichlet.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p)
