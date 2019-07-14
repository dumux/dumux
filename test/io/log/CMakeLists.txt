set(CMAKE_BUILD_TYPE Debug)

dune_add_test(SOURCES test_log.cc
              LABELS unit io)
