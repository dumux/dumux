add_input_file_links()
dune_symlink_to_source_files(FILES grids)

# build target for the energy kinetic test problem, two energy balance equations
dune_add_test(COMPILE_ONLY # since it currently fails miserably with very different results on different machines
              NAME test_mpnc_thermalnonequil_box
              SOURCES main.cc
              LABELS mpnc
              COMPILE_DEFINITIONS TYPETAG=CombustionOneComponentBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_thermalnonequil_box-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_thermalnonequil_box-00073.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_thermalnonequil_box params.input -Problem.Name test_mpnc_thermalnonequil_box")

#install sources
install(FILES
combustionfluidsystem.hh
combustionlocalresidual.hh
problem.hh
spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/mpnc/thermalnonequilibrium)
