add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dune_add_test(NAME test_dec1p
               SOURCES test_1p.cc
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p-00000.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_dec1p -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_1p.input")

dune_add_test(NAME test_diffusion
              SOURCES test_diffusion.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/diffusion-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/mimeticdiffusion-00000.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_diffusion")

dune_add_test(NAME test_diffusion3d
              SOURCES test_diffusion3d.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_diffusion3d_fv-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion3d_fv-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_diffusion3d_fvmpfal-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion3d_fvmpfal-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_diffusion3d_mimetic-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion3d_mimetic-00000.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_diffusion3d -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_diffusion3d_reference.input")

#install sources
install(FILES
resultevaluation3d.hh
resultevaluation.hh
test_1p.cc
test_1pproblem.hh
test_1pspatialparams.hh
test_diffusion3d.cc
test_diffusion.cc
test_diffusionproblem3d.hh
test_diffusionproblem.hh
test_diffusionspatialparams3d.hh
test_diffusionspatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/1p/sequential)
