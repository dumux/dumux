add_input_file_links()
dune_symlink_to_source_files(FILES generate_network_tubes.py calculate_analytical_solution.py)
dumux_add_test(NAME validate_preferential_path
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS CIRCLETHROAT=1
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=0
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME validate_preferential_path_kilian
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS CIRCLETHROAT=1
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")
