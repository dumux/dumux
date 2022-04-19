dune_symlink_to_source_files(FILES "convergencetest.py")
dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_ff_navierstokesnc_analytical
              SOURCES main.cc
              LABELS freeflow  navierstokes navierstokesnc analytical
              TIMEOUT 3600
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS test_ff_navierstokesnc_analytical params.input
              -Grid.LowerLeft "0.0 0.0"
              -Grid.UpperRight "1.0 1.0"
              -Grid.Cells "25 25"
              -Problem.Name test_ff_navierstokesnc_analytical
              -Problem.PrintConvergenceTestFile true)
