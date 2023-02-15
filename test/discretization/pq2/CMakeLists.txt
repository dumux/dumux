dumux_add_test(NAME test_pq2fvgeometry
               SOURCES test_pq2fvgeometry.cc
               COMPILE_DEFINITIONS ENABLE_CACHING=true
               LABELS unit discretization)

set(CMAKE_BUILD_TYPE Debug)