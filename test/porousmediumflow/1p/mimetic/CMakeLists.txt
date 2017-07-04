add_input_file_links()

# isothermal tests
add_dumux_test(test_mimetic1p test_mimetic1p test_mimetic1p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_mimetic1p")

add_dumux_test(test_1pgermanbasin test_1pgermanbasin test_1pgermanbasin.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pgermanbasin)

add_dumux_test(test_1pmimeticanisotropic test_1pmimeticanisotropic test_1pmimeticanisotropic.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pmimeticanisotropic)

#install sources
install(FILES
1ptestproblem.hh
1ptestspatialparams.hh
tubesspatialparams.hh
test_mimetic1p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
