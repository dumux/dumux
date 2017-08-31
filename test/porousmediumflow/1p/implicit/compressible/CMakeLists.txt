dune_symlink_to_source_files(FILES "test_1p.input")

# compressible
dune_add_test(NAME test_1p_compressible_tpfa
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleTpfa
              CMD_ARGS test_1p.input -Problem.Name 1p_tpfa)

dune_add_test(NAME test_1p_compressible_box
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleBox
              CMD_ARGS test_1p.input -Problem.Name 1p_box)

set(CMAKE_BUILD_TYPE Release)
