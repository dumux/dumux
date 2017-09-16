dune_symlink_to_source_files(FILES "test_1p.input" "test_1p_stationary.input")

# compressible
dune_add_test(NAME test_1p_compressible_tpfa
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleTpfa
              CMD_ARGS test_1p.input -Problem.Name 1p_tpfa)

dune_add_test(NAME test_1p_compressible_mpfa
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleMpfa
              CMD_ARGS test_1p.input -Problem.Name 1p_mpfa)

dune_add_test(NAME test_1p_compressible_box
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleBox
              CMD_ARGS test_1p.input -Problem.Name 1p_box)

# stationary
dune_add_test(NAME test_1p_compressible_stationary_tpfa
              SOURCES test_1p_stationary.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleTpfa
              CMD_ARGS test_1p_stationary.input -Problem.Name 1p_stationary_tpfa)

dune_add_test(NAME test_1p_compressible_stationary_mpfa
              SOURCES test_1p_stationary.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleMpfa
              CMD_ARGS test_1p_stationary.input -Problem.Name 1p_stationary_mpfa)

dune_add_test(NAME test_1p_compressible_stationary_box
              SOURCES test_1p_stationary.cc
              COMPILE_DEFINITIONS TYPETAG=OnePCompressibleBox
              CMD_ARGS test_1p_stationary.input -Problem.Name 1p_stationary_box)

set(CMAKE_BUILD_TYPE Release)
