add_input_file_links()

dune_add_test(NAME test_stokes1pdarcy1pconv
              SOURCES test_stokes1pdarcy1pconv.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS DARCYGRIDTYPE=Dune::YaspGrid<2>
              COMMAND ./test_stokes1pdarcy1pconv
              CMD_ARGS "test_stokes1pdarcy1pconv.input")
