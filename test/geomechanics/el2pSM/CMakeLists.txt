add_input_file_links()

# for the test only simulate one week
add_dumux_test(test_el2pdecoupledSM test_el2pdecoupledSM test_el2pdecoupledSM.cc echo)

# headers for installation and headercheck
install(FILES
        decoupled2pproblemLu.hh
        decoupled2pspatialparamsLu.hh
        decoupledelastic2pproblemLu.hh
        decoupledelasticproblemLu.hh
        test_el2pdecoupledSM.cc
        test_el2pdecoupledSM.input
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/geomechanics/el2pSM)
