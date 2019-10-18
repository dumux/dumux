add_input_file_links()

dune_symlink_to_source_files(FILES grids)

# dim < dimWorld tests with Dune::Foamgrid<2,3>
dumux_add_test(NAME test_1p_fracture2d3d_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FractureBox
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_fracture2d3d_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_fracture2d3d_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_fracture2d3d_box params.input -Problem.Name test_1p_fracture2d3d_box")

dumux_add_test(NAME test_1p_fracture2d3d_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FractureCCTpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_fracture2d3d_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_fracture2d3d_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_fracture2d3d_tpfa params.input -Problem.Name test_1p_fracture2d3d_tpfa")

dumux_add_test(NAME test_1p_fracture2d3d_mpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FractureCCMpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_fracture2d3d_cc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_1p_fracture2d3d_mpfa-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_fracture2d3d_mpfa params.input -Problem.Name test_1p_fracture2d3d_mpfa")
