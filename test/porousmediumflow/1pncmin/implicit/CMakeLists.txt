# isothermal tests

add_input_file_links()

# add_executable(test_box1pncmin test_box1pncmin.cc)
# add_executable(test_cc1pncmin test_cc1pncmin.cc)

add_dumux_test(test_box1pncmin test_box1pncmin test_box1pncmin.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/box1pncmin-00042.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1pncmin-00042.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pncmin"
                 --zeroThreshold {"temperature":5e-6})

add_dumux_test(test_cc1pncmin test_cc1pncmin test_cc1pncmin.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/cc1pncmin-00042.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1pncmin-00042.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pncmin"
                 --zeroThreshold {"temperature":5e-6})
#
# #install sources
install(FILES
thermochempatialparams.hh
thermochemproblem.hh
test_box1pncmin.cc
test_cc1pncmin.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1pncmin)
