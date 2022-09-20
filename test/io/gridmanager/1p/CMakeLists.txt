add_input_file_links(params.input)

dune_add_test(NAME test_1p_landfillgridmanager_box
              LABELS porousmediumflow 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=LandfillGridmanagerTestBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py)

dune_add_test(NAME test_1p_landfillgridmanager_tpfa
              LABELS porousmediumflow 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=LandfillGridmanagerTestCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py)
