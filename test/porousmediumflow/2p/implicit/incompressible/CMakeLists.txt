dune_symlink_to_source_files(FILES "test_2p.input")

# using tpfa
dune_add_test(NAME test_2p_incompressible_tpfa
              SOURCES test_2p.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleTpfa
              CMD_ARGS test_2p.input -Problem.Name 2p_tpfa)

# using box
dune_add_test(NAME test_2p_incompressible_box
              SOURCES test_2p.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleBox
              CMD_ARGS test_2p.input -Problem.Name 2p_box)

set(CMAKE_BUILD_TYPE Release)
