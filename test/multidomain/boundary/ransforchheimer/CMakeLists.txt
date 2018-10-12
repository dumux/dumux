dune_add_test(NAME test_rans1pdarcy1pforchheimer
              SOURCES test_rans1pdarcy1pforchheimer.cc
              COMPILE_DEFINITIONS FORCHHEIMER=1
              CMAKE_GUARD HAVE_UMFPACK
              CMD_ARGS test_rans1pdarcy1pforchheimer.input)

dune_symlink_to_source_files(FILES "test_rans1pdarcy1pforchheimer.input")
