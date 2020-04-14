dune_symlink_to_source_files(FILES config.ini coordinates.txt mydata.xml numbers.txt)
dumux_add_test(SOURCES test_io_data_input.cc LABELS unit io)
