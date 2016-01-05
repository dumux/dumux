add_input_file_links()

add_dumux_test(test_cc1p_pointsources test_cc1p_pointsources test_cc1p_pointsources.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1psingularitycc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1psingularitycc-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p_pointsources")

add_dumux_test(test_box1p_pointsources test_box1p_pointsources test_box1p_pointsources.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1psingularitybox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1psingularitybox-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p_pointsources")


#install sources
install(FILES
1psingularityproblem.hh
1psingularityspatialparams.hh
test_cc1p_pointsources.cc
test_box1p_pointsources.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/1p/implicit/pointsources)
