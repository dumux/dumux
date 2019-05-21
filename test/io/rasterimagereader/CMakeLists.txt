dune_symlink_to_source_files(FILES blackwhite_j.pbm blackwhite_binary_j.pbm blackwhite_j.txt grayscale_j.pgm grayscale_binary_j.pgm grayscale_j.txt)

dune_add_test(NAME test_io_rasterimagereader
              SOURCES test_rasterimagereader.cc
              LABELS unit io)
