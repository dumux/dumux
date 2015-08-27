add_input_file_links()

# isothermal tests
add_dumux_test(test_box1p test_box1p test_box1p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p")

add_dumux_test(test_box1pwithamg test_box1pwithamg test_box1pwithamg.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestboxwithamg-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pwithamg")

add_dumux_test(test_cc1p test_cc1p test_cc1p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p")

add_dumux_test(test_cc1pwithamg test_cc1pwithamg test_cc1pwithamg.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestccwithamg-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pwithamg")

# non-isothermal tests
add_dumux_test(test_box1pniconduction test_box1pniconduction test_box1pniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1pniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pniconduction"
                 --zeroThreshold {"velocity":1e-8})

add_dumux_test(test_box1pniconvection test_box1pniconvection test_box1pniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1pniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pniconvection"
                 --zeroThreshold {"velocity":1e-16})

add_dumux_test(test_cc1pniconduction test_cc1pniconduction test_cc1pniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pniccconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1pniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pniconduction"
                 --zeroThreshold {"velocity":1e-8})

add_dumux_test(test_cc1pniconvection test_cc1pniconvection test_cc1pniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1pniccconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1pniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pniconvection")

# dim < dimWorld tests
add_dumux_test(test_box1p1d3d test_box1p1d3d test_box1p1d3d.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox1d3d-reference.vtp
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox1d3d-00002.vtp
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p1d3d")
add_dune_alberta_flags(test_box1p1d3d GRIDDIM 1 WORLDDIM 3)

add_dumux_test(test_box1p2d3d test_box1p2d3d test_box1p2d3d.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox2d3d-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox2d3d-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2d3d")
add_dune_alberta_flags(test_box1p2d3d GRIDDIM 2 WORLDDIM 3)

add_dumux_test(test_cc1p1d3d test_cc1p1d3d test_cc1p1d3d.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc1d3d-reference.vtp
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc1d3d-00002.vtp
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p1d3d")
add_dune_alberta_flags(test_cc1p1d3d GRIDDIM 1 WORLDDIM 3)

add_dumux_test(test_cc1p2d3d test_cc1p2d3d test_cc1p2d3d.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc2d3d-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc2d3d-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2d3d")
add_dune_alberta_flags(test_cc1p2d3d GRIDDIM 2 WORLDDIM 3)


#install sources
install(FILES
1pniconductionproblem.hh
1pniconvectionproblem.hh
1pnispatialparams.hh
1ptestproblem.hh
1ptestspatialparams.hh
test_box1p1d3d.cc
test_box1p2d3d.cc
test_box1p.cc
test_box1pniconduction.cc
test_box1pniconvection.cc
test_box1pwithamg.cc
test_cc1p1d3d.cc
test_cc1p2d3d.cc
test_cc1p.cc
test_cc1pniconduction.cc
test_cc1pniconvection.cc
test_cc1pwithamg.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
