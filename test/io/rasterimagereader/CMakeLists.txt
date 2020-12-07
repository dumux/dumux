dune_symlink_to_source_files(FILES blackwhite_j.pbm blackwhite_binary_j.pbm blackwhite_j.txt
                                   grayscale_j.pgm grayscale_binary_j.pgm grayscale_j.txt
                                   blackwhite_dim_firstline.pbm blackwhite_dim_firstline_binary.pbm
                                   blackwhite_fail.pbm grayscale_fail_binary.pgm)

dumux_add_test(NAME test_io_rasterimagereader
              SOURCES test_rasterimagereader.cc
              LABELS unit io)
