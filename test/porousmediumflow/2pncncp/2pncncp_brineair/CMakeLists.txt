add_input_file_links()

dumux_add_test(NAME test_twopnc_ncp_brineair_box
              LABELS porousmediumflow
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPNCSaltBox)

dumux_add_test(NAME test_twopnc_ncp_brineair_tpfa
              LABELS porousmediumflow
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPNCSaltCC)
