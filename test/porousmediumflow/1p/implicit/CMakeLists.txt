add_subdirectory(pointsources)

add_input_file_links()
add_gstat_file_links()

dune_symlink_to_source_files(FILES "tubesconvergencetest.py")

add_executable(test_mpfa1p EXCLUDE_FROM_ALL test_mpfa1p.cc)

# isothermal tests
add_dumux_test(test_box1p test_box1p test_box1p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p")

add_dumux_test(test_box1pwithamg test_box1pwithamg test_box1pwithamg.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestboxwithamg-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pwithamg")

add_dumux_test(test_cc1p test_cc1p test_cc1p.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p")

add_dumux_test(test_cc1pwithamg test_cc1pwithamg test_cc1pwithamg.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestccwithamg-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pwithamg")

add_dumux_test(test_cc1pwithgstat test_cc1pwithgstat test_cc1pwithgstat.cc
               ./test_cc1pwithgstat)

# non-isothermal tests
add_dumux_test(test_box1pniconduction test_box1pniconduction test_box1pniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1pniconduction-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pniconduction"
                 --zeroThreshold {"velocity":1e-8})

add_dumux_test(test_box1pniconvection test_box1pniconvection test_box1pniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1pniconvection-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pniconvection"
                 --zeroThreshold {"velocity":1e-15})

add_dumux_test(test_cc1pniconduction test_cc1pniconduction test_cc1pniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pniccconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1pniconduction-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pniconduction"
                 --zeroThreshold {"velocity":1e-8})

add_dumux_test(test_cc1pniconvection test_cc1pniconvection test_cc1pniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pniccconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1pniconvection-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pniconvection")

# dim < dimWorld tests with foamgrid
add_dumux_test(test_box1pnetwork1d3d test_box1pnetwork1d3d test_box1pnetwork1d3d.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pboxnetwork1d3d-reference.vtp
                         ${CMAKE_CURRENT_BINARY_DIR}/1pboxnetwork1d3d-00001.vtp
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pnetwork1d3d")

add_dumux_test(test_box1pfracture2d3d test_box1pfracture2d3d test_box1pfracture2d3d.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pboxfracture2d3d-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1pboxfracture2d3d-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pfracture2d3d")

dune_add_test(SOURCES test_cc1pnetwork1d3d.cc
              NAME test_cc1pnetwork1d3d
              COMMAND ./tubesconvergencetest.py test_cc1pnetwork1d3d)

add_dumux_test(test_cc1pfracture2d3d test_cc1pfracture2d3d test_cc1pfracture2d3d.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pccfracture2d3d-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1pccfracture2d3d-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pfracture2d3d")


#install sources
install(FILES
1pniconductionproblem.hh
1pniconvectionproblem.hh
1pnispatialparams.hh
1ptestproblem.hh
1ptestspatialparams.hh
1pboxfracture2d3d.cc
1pboxnetwork1d3d.cc
test_box1p.cc
test_box1pniconduction.cc
test_box1pniconvection.cc
test_box1pwithamg.cc
test_cc1pfracture2d3d.cc
test_cc1pnetwork1d3d.cc
test_cc1p.cc
test_cc1pniconduction.cc
test_cc1pniconvection.cc
test_cc1pwithamg.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
