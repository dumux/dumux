add_input_file_links()
dune_symlink_to_source_files(FILES dt_vs_efficiency_accruacy.py run_and_plot_convergence_drainage.py)
dune_symlink_to_source_files(FILES run_and_plot_convergence_imbibition.py extractresults.py plot_constitutive_laws.py)

dumux_add_test(NAME test_pnm_2p_1d_oilwater_drainage
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=1
               COMPILE_DEFINITIONS CIRCLETHROAT=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME test_pnm_2p_1d_oilwater_drainage_other_implicit
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=1
               COMPILE_DEFINITIONS CIRCLETHROAT=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")



dumux_add_test(NAME test_pnm_2p_1d_oilwater_imbibition
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=0
               COMPILE_DEFINITIONS ONEDDRAINGE=1
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")


dumux_add_test(NAME test_pnm_2p_1d_oilwater_imbibition_other_implicit
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=0
               COMPILE_DEFINITIONS ONEDDRAINGE=1
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMPILE_DEFINITIONS NOREGULARIZTAION=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME plot_constitutive_relations
               SOURCES main_plot_constitutive_laws.cc
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMPILE_DEFINITIONS DRAINAGE=0
               COMPILE_DEFINITIONS ONEDDRAINGE=1
               COMPILE_DEFINITIONS ENABLECACHE=0
               COMPILE_DEFINISIONS REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS PLOTCONSTITUTIVERELATIONS=1)
