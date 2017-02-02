add_input_file_links()


add_dumux_test(test_boxrichardsnc test_boxrichardsnc test_boxrichardsnc.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardswelltracerbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/richardswelltracerbox-00032.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsnc")

add_dumux_test(test_ccrichardsnc test_ccrichardsnc test_ccrichardsnc.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardswelltracercc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/richardswelltracercc-00031.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsnc")

#install sources
install(FILES
richardswelltracerproblem.hh
richardswelltracerspatialparams.hh
test_boxrichardsnc.cc
test_ccrichardsnc.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/richardsnc/implicit)
