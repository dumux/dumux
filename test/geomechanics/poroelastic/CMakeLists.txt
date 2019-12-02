dune_symlink_to_source_files(FILES "params.input")

# using box and numeric differentiation
dumux_add_test(NAME test_poroelastic_box
              LABELS geomechanics poroelastic
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_poroelastic_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_box params.input -Problem.Name test_poroelastic_box")
