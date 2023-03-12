add_input_file_links()
dune_symlink_to_source_files(FILES grids plot_curve.gp)

dune_add_test(NAME pc_sw_static
              SOURCES main_static.cc)

dune_add_test(NAME pc_sw_dynamic
              SOURCES main_dynamic.cc)
