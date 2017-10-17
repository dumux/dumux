# isothermal tests

add_input_file_links()

add_executable(test_box1pncmin test_box1pncmin.cc)
add_executable(test_cc1pncmin test_cc1pncmin.cc)

# add_executable_all(test_box1p2c test_box1p2c.cc)
#
# add_executable_all(test_cc1p2c test_cc1p2c.cc)

# dune_symlink_to_source_files(FILES test_2pnc.input)
# # dune_symlink_to_source_files(FILES test_box2pnc.input)
# dune_symlink_to_source_files(FILES test_box1p2c.input)
# dune_symlink_to_source_files(FILES test_cc1p2c.input)

# set(BUILD_TYPE Debug)

#install sources
# install(FILES
# thermochemproblem.hh
# thermochemspatialparams.hh
# test_box2pnc.cc
# test_cc2pnc.cc
# DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-devel/appl/thermochemistry)


# # isothermal tests
# add_dumux_test(test_box2pnc test_box2pnc test_box2pnc.cc
#                python ${dumux_INCLUDE_DIRS}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pncbox-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/fuelcell-box-00016.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pnc -ParameterFile test_2pnc.input -Problem.Name fuelcell-box")
#
# add_dumux_test(test_cc2pnc test_cc2pnc test_cc2pnc.cc
#                python ${dumux_INCLUDE_DIRS}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pnccc-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/fuelcell-cc-00016.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pnc -ParameterFile test_2pnc.input -Problem.Name fuelcell-cc")
#
# dune_symlink_to_source_files(FILES test_2pnc.input)
#
# #install sources
# install(FILES
# thermochemproblem.hh
# thermochemspatialparams.hh
# test_box2pnc.cc
# test_cc2pnc.cc
# DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pnc)
