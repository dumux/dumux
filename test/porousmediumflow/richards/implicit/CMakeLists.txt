add_input_file_links()

if(MPI_FOUND)
  add_dumux_test(test_boxrichards test_boxrichards test_boxrichards.cc
                 python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                   --script fuzzy
                   --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference-parallel.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-richardslensbox-00008.vtu
                   --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards")
else()
# isothermal tests
  add_dumux_test(test_boxrichards test_boxrichards test_boxrichards.cc
                 python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                   --script fuzzy
                   --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/richardslensbox-00008.vtu
                   --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards")
endif()

add_dumux_test(test_ccrichards test_ccrichards test_ccrichards.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/richardslenscc-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichards")

# comparison to analytical solution - only with cc
add_dumux_test(test_ccrichardsanalytical test_ccrichardsanalytical test_ccrichardsanalytical.cc
                 python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardsanalyticalcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/richardsanalyticalcc-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsanalytical")

# non-isothermal tests
add_dumux_test(test_boxrichardsniconvection test_boxrichardsniconvection test_boxrichardsniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection"
                 --zeroThreshold {"velocity":1e-16})

add_dumux_test(test_ccrichardsniconvection test_ccrichardsniconvection test_ccrichardsniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectioncc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconvection")

add_dumux_test(test_boxrichardsniconduction test_boxrichardsniconduction test_boxrichardsniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconductionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction"
                 --zeroThreshold {"velocity":1e-8})

add_dumux_test(test_ccrichardsniconduction test_ccrichardsniconduction test_ccrichardsniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconductioncc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction"
                 --zeroThreshold {"velocity":1e-8})

#install sources
install(FILES
richardsanalyticalproblem.hh
richardsanalyticalspatialparams.hh
richardslensproblem.hh
richardslensspatialparams.hh
richardsniconductionproblem.hh
richardsniconvectionproblem.hh
richardsnispatialparams.hh
test_boxrichards.cc
test_boxrichardsniconduction.cc
test_boxrichardsniconvection.cc
test_ccrichardsanalytical.cc
test_ccrichards.cc
test_ccrichardsniconduction.cc
test_ccrichardsniconvection.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/richards)
