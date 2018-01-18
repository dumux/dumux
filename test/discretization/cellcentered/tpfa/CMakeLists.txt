dune_add_test(NAME test_tpfafvgeometry
              SOURCES test_tpfafvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false)

dune_add_test(NAME test_tpfafvgeometry_caching
              SOURCES test_tpfafvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true)

dune_add_test(NAME test_tpfafvgeometry_nonconforming
              SOURCES test_tpfafvgeometry_nonconforming.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false
              CMAKE_GUARD dune-alugrid_FOUND)

dune_add_test(NAME test_cachedtpfafvgeometry_nonconforming
              SOURCES test_tpfafvgeometry_nonconforming.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true
              CMAKE_GUARD dune-alugrid_FOUND)
