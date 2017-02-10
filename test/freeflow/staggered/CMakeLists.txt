add_input_file_links()

add_dumux_test(test_liddrivencavity test_closedsystem test_closedsystem.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/liddrivencavity-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_liddrivencavity-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_liddrivencavity.input")

add_dumux_test(test_hydrostaticpressure test_closedsystem test_closedsystem.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/stokeshydrostaticpressure-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_hydrostaticpressure-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_hydrostaticpressure.input"
                 --zeroThreshold {"velocity_Constant \(m/s\)":1e-16})

add_dumux_test(test_channel_stokes test_channel_stokes test_channel.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/channel-stokes.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokes-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_stokes")

add_dumux_test(test_channel_navierstokes test_channel_navierstokes test_channel.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/channel-navierstokes-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_channel_navierstokes-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_navierstokes")
target_compile_definitions(test_channel_navierstokes PUBLIC "ENABLE_NAVIERSTOKES=1")

add_dumux_test(test_donea test_donea test_donea.cc
      python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
        --script fuzzy
        --files ${CMAKE_SOURCE_DIR}/test/references/stokes-donea-reference.vtu
                ${CMAKE_CURRENT_BINARY_DIR}/test_donea-00002.vtu
        --command "${CMAKE_CURRENT_BINARY_DIR}/test_donea")

add_dumux_test(test_kovasznay test_kovasznay test_kovasznay.cc
      python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
        --script fuzzy
        --files ${CMAKE_SOURCE_DIR}/test/references/test_kovasznay-reference.vtu
                ${CMAKE_CURRENT_BINARY_DIR}/test_kovasznay-00002.vtu
        --command "${CMAKE_CURRENT_BINARY_DIR}/test_kovasznay")

#install sources
install(FILES
channeltestproblem.hh
test_channel.cc
closedsystemtestproblem.hh
test_closedsystem.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/staggered)
