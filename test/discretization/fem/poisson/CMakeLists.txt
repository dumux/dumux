dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_fem_poisson2d
              SOURCES test_fem_poisson2d.cc
              COMPILE_DEFINITIONS TENSORTYPE=Scalar
              LABELS unit discretization fem)
set(CMAKE_BUILD_TYPE Release)
