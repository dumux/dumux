dune_symlink_to_source_files(FILES "test_1p.input")

# using tpfa
dune_add_test(NAME test_1p_incompressible_tpfa
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa
              CMD_ARGS test_1p.input)

# using mpfa
dune_add_test(NAME test_1p_incompressible_mpfa
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleMpfa
              CMD_ARGS test_1p.input)

# using box
dune_add_test(NAME test_1p_incompressible_box
              SOURCES test_1p.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox
              CMD_ARGS test_1p.input)

set(CMAKE_BUILD_TYPE Release)
