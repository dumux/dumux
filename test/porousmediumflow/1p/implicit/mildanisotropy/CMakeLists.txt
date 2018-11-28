dune_symlink_to_source_files(FILES "params.input")

# using tpfa
dune_add_test(NAME test_1p_mildanisotropy_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MildAnisotropyTpfa
              COMMAND ./test_1test_1p_mildanisotropy_tpfa
              CMD_ARGS  "params.input")

# using mpfa
dune_add_test(NAME test_1p_mildanisotropy_mpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MildAnisotropyMpfa
              COMMAND ./test_1test_1p_mildanisotropy_mpfa
              CMD_ARGS  "params.input")

# using box
dune_add_test(NAME test_1p_mildanisotropy_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MildAnisotropyBox
              COMMAND ./test_1test_1p_mildanisotropy_box
              CMD_ARGS  "params.input")

set(CMAKE_BUILD_TYPE Release)
