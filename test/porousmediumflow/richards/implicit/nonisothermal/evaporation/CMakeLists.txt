add_input_file_links(FILES params.input)

# non-isothermal tests
dune_add_test(NAME test_richardsni_evaporation_tpfa
              LABELS porousmediumflow richards richardsni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsNIEvaporationCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richardsni_evaporation_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_evaporation_tpfa-00043.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_evaporation_tpfa params.input -Problem.Name test_richardsni_evaporation_tpfa")

dune_add_test(NAME test_richardsni_evaporation_box
              LABELS porousmediumflow richards richardsni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsNIEvaporationBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richardsni_evaporation_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_evaporation_box-00043.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_evaporation_box params.input -Problem.Name test_richardsni_evaporation_box")
