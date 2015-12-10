add_input_file_links()

# isothermal tests
add_dumux_test(test_box3p test_box3p test_box3p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/infiltration3pbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pbox-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3p")

add_dumux_test(test_cc3p test_cc3p test_cc3p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/infiltration3pcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pcc-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3p")

# non-isothermal tests
add_dumux_test(test_box3pniconvection test_box3pniconvection test_box3pniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/3pniconvectionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconvection"
                 --zeroThreshold {"velocityW":1e-16})

add_dumux_test(test_cc3pniconvection test_cc3pniconvection test_cc3pniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/3pniconvectioncc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconvection"
                 --zeroThreshold {"velocityW":1e-16})

add_dumux_test(test_box3pniconduction test_box3pniconduction test_box3pniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/3pniconductionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconduction"
                 --zeroThreshold {"velocityW":1e-8})

add_dumux_test(test_cc3pniconduction test_cc3pniconduction test_cc3pniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/3pniconductioncc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconduction"
                 --zeroThreshold {"velocityW":1e-8})

#install sources
install(FILES
3pniconductionproblem.hh
3pniconvectionproblem.hh
3pnispatialparams.hh
infiltration3pproblem.hh
infiltration3pspatialparams.hh
test_box3p.cc
test_box3pniconduction.cc
test_box3pniconvection.cc
test_cc3p.cc
test_cc3pniconduction.cc
test_cc3pniconvection.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/3p)
