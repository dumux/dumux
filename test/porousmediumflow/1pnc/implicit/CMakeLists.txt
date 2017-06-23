add_input_file_links()

# isothermal tests
add_executable(test test_box1p2c.cc)

# add_executable(test_cc1pnc test_cc1pnc.cc)




# add_dumux_test(test_box1p2c test_box1p2c test_box1p2c.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/outflowbox-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/outflowbox-00009.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2c"
#                  --zeroThreshold {"velocity":5e-16})
#
# add_dumux_test(test_cc1p2c test_cc1p2c test_cc1p2c.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/1p2ctestcc-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/1p2ctestcc-00009.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2c"
#                  --zeroThreshold {"velocity":5e-16})
#
# #install sources
# install(FILES
#
# 1p2cnispatialparams.hh
# 1p2coutflowproblem.hh
# 1p2coutflowspatialparams.hh
# test_box1p2c.cc
# test_box1p2cniconduction.cc
# test_box1p2cniconvection.cc
# test_cc1p2c.cc
# test_cc1p2cniconduction.cc
# test_cc1p2cniconvection.cc
# DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1pnc)
