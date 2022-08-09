dumux_add_test(NAME test_cvfefvgeometry
              SOURCES test_cvfefvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true
              LABELS unit discretization)
