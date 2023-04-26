dune_symlink_to_source_files(FILES "convergencetest.py")
dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_ff_navierstokesnc_sincos_oldstaggered
              SOURCES main.cc
              LABELS freeflow  navierstokes navierstokesnc
              TIMEOUT 3600
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS LINEARSOLVER=UMFPackBackend
              COMMAND ./convergencetest.py
              CMD_ARGS test_ff_navierstokesnc_sincos params.input
              -Grid.UpperRight "6.28 6.28"
              -Grid.Cells "50 50"
              -Problem.Name test_ff_navierstokesnc_sincos
              -Problem.PrintConvergenceTestFile true)
