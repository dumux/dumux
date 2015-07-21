add_subdirectory("2p")

add_dumux_test(test_effectivediffusivitymillingtonquirk test_effectivediffusivitymillingtonquirk test_effectivediffusivitymillingtonquirk.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script exact
                 --files ${CMAKE_SOURCE_DIR}/test/references/effectivediffusivitymillingtonquirk-reference.dat
                         ${CMAKE_CURRENT_BINARY_DIR}/millingtonquirk_d_eff.dat
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_effectivediffusivitymillingtonquirk")

add_dumux_test(test_effectivediffusivityconstanttau test_effectivediffusivityconstanttau test_effectivediffusivityconstanttau.cc
              python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script exact
                 --files ${CMAKE_SOURCE_DIR}/test/references/effectivediffusivityconstanttau-reference.dat
                         ${CMAKE_CURRENT_BINARY_DIR}/constanttau_d_eff.dat
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_effectivediffusivityconstanttau")
