file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
dumux_add_test(SOURCES test_gnuplotinterface.cc
              LABELS unit io)
