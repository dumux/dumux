dune_add_test(NAME test_tpfafvgeometry
              SOURCES test_tpfafvgeometry.cc
              COMMAND ./test_tpfafvgeometry)

dune_add_test(NAME test_tpfafvgeometry_nonconforming
              SOURCES test_tpfafvgeometry_nonconforming.cc
              COMPILE_DEFINITIONS TYPETAG=TestFVGeometryNonConforming
              COMMAND ./test_tpfafvgeometry_nonconforming)

dune_add_test(NAME test_cachedtpfafvgeometry_nonconforming
              SOURCES test_tpfafvgeometry_nonconforming.cc
              COMPILE_DEFINITIONS TYPETAG=TestCachedFVGeometryNonConforming
              COMMAND ./test_cachedtpfafvgeometry_nonconforming)
