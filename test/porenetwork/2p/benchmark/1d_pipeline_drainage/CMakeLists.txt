add_input_file_links()
dune_symlink_to_source_files(FILES check_invading_time.py plot_PcS_and_derivatives.py)

dumux_add_test(NAME test_pnm_2p_reg
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=1
               COMPILE_DEFINITIONS CIRCLETHROAT=1
               COMPILE_DEFINITIONS STOPDRAINAGEAFTERFIRSTTHROAT=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME test_pnm_2p_no_reg
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=1
               COMPILE_DEFINITIONS CIRCLETHROAT=1
               COMPILE_DEFINITIONS STOPDRAINAGEAFTERFIRSTTHROAT=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")
