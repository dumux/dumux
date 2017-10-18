add_input_file_links()

# add_executable(test_box1pnc test_box1pnc.cc)
# add_executable(test_cc1pnc test_cc1pnc.cc)

add_dumux_test(test_box1pnc test_box1pnc test_box1pnc.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pnctestbox-00009.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1pnctestbox-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pnc"
                 --zeroThreshold {"pressure":5e-6})

add_dumux_test(test_cc1pnc test_cc1pnc test_cc1pnc.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pnctestcc-00009.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1pnctestcc-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pnc"
                 --zeroThreshold {"pressure":5e-6})
#
# #install sources
install(FILES
1pncspatialparams.hh
1pnctestproblem.hh
test_box1pnc.cc
test_cc1pnc.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1pnc)
