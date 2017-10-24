add_subdirectory("convergence")
add_subdirectory("DMP")

add_input_file_links()

dune_symlink_to_source_files(FILES runAllTestsWeakConv.sh)
dune_symlink_to_source_files(FILES runTestWeakConv.sh)
dune_symlink_to_source_files(FILES runTestWeakConvMPFAL.sh)

# isothermal tests
add_dumux_test(test_mimetic1p test_mimetic1p test_mimetic1p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_mimetic1p")

add_dumux_test(test_1pgermanbasin test_1pgermanbasin test_1pgermanbasin.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pgermanbasin)

add_dumux_test(test_1plinear test_1plinear test_1plinear.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1plinear)

add_dumux_test(test_mimetic1pniconduction test_mimetic1pniconduction test_mimetic1pniconduction.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_mimetic1pniconduction)

add_dumux_test(test_1pweakconv test_1pweakconv test_1pweakconv.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pweakconv)

add_dumux_test(test_1pweakconv_mpfal test_1pweakconv_mpfal test_1pweakconv_mpfal.cc
                ${CMAKE_CURRENT_BINARY_DIR}/test_1pweakconv_mpfal)

#install sources
install(FILES
1ptestproblem.hh
1ptestspatialparams.hh
tubesspatialparams.hh
test_mimetic1p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
