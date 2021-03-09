dune_symlink_to_source_files(FILES params.input test_1p2c_saltwaterintrusion_box_dt-reference.dat)

# salt-water intrusion test using the box model
dumux_add_test(NAME test_1p2c_saltwaterintrusion_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=SaltWaterIntrusionTest
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_saltwaterintrusion_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_saltwaterintrusion_box-00033.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_saltwaterintrusion_box params.input -Problem.Name test_1p2c_saltwaterintrusion_box -TimeLoop.TimeStepSizeFile test_1p2c_saltwaterintrusion_box_dt-reference.dat")
