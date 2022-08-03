dumux_add_test(NAME test_diamondgridgeometry
               SOURCES test_diamondgridgeometry.cc
               COMPILE_DEFINITIONS ENABLE_CACHING=true
               LABELS unit discretization)
