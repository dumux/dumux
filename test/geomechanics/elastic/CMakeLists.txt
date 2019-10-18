dune_symlink_to_source_files(FILES "params.input")

# using box
dumux_add_test(NAME test_elastic_box
              LABELS geomechanics elastic
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_elastic_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_elastic_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_elastic_box params.input -Problem.Name test_elastic_box")
