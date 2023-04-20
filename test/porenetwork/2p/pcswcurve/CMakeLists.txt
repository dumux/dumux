add_input_file_links()
dune_symlink_to_source_files(FILES grids run_and_plot.py)

dune_add_test(NAME pc_sw_static
              SOURCES main_static.cc)

dune_add_test(NAME pc_sw_dynamic
              SOURCES main_dynamic.cc
              COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
              COMPILE_DEFINITIONS CIRCLETHROAT=0
              COMIPLE_DEFINITIONS DRAINAGE=1
              COMPILE_DEFINITIONS NOREGULARIZTAION=0)

dune_add_test(NAME pc_sw_dynamic_other_implict
              SOURCES main_dynamic.cc
              COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
              COMPILE_DEFINITIONS CIRCLETHROAT=0
              COMPILE_DEFINITIONS DRAINAGE=1
              COMPILE_DEFINITIONS NOREGULARIZTAION=1)
