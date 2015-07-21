if(MPI_CXX_FOUND)
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
endif(MPI_CXX_FOUND)

add_dumux_test(test_ccrichards test_ccrichards test_ccrichards.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/richardslenscc-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichards")

# non-isothermal tests
add_dumux_test(test_boxrichardsniconvection test_boxrichardsniconvection test_boxrichardsniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection")

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
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction")

add_dumux_test(test_ccrichardsniconduction test_ccrichardsniconduction test_ccrichardsniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconductioncc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction")