add_input_file_links()
dune_symlink_to_source_files(FILES 4_homo_pores.dgf 8-pores-network.dgf run_and_plot.py iterations_and_timesteps.py)

dumux_add_test(NAME pnm_3d_comparison_manchester
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=0
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME pnm_3d_comparison_no_regularization_manchester
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=1
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")
