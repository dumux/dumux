dune_symlink_to_source_files(FILES bspline_plot.py)

dumux_add_test(SOURCES test_function_l2norm.cc
              CMAKE_GUARD dune-functions_FOUND
              LABELS unit)

dumux_add_test(SOURCES test_bspline_function.cc
               CMAKE_GUARD dune-functions_FOUND
               LABELS unit)
