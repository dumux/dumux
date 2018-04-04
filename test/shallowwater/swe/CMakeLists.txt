dune_add_test(NAME test_swe
              SOURCES test_swe.cc
              COMPILE_DEFINITIONS TYPETAG=SweTestTypeTag
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_swe)
