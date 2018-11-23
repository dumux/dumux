add_input_file_links()

dune_add_test(NAME test_1pncminni_charge_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ChargeBox
              CMD_ARGS params.input)

#install sources
install(FILES
problem.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/1pncmin/charge)
