add_input_file_links()

# salt-water intrusion test using the box model
dune_add_test(NAME test_1p2c_saltwaterintrusion_box
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=SaltWaterIntrusionTest
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_saltwaterintrusion_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_saltwaterintrusion_box-00035.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_saltwaterintrusion_box params.input -Problem.Name test_1p2c_saltwaterintrusion_box")
