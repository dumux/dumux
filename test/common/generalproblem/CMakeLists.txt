add_input_file_links()

add_dumux_test(test_general_box test_general test_generalproblem2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/generallens_box-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/generallens_box-00003.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_general -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p_reference.input -ModelType box")

add_dumux_test(test_general_cc test_general test_generalproblem2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/generallens_cc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/generallens_cc-00003.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_general -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p_reference.input -ModelType cc")

add_dumux_test(test_general_dec test_general test_generalproblem2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/generallens_decoupled-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/generallens_decoupled-00003.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_general -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_generalproblem2p_reference.input -ModelType sequential")

#install sources
install(FILES
generallensproblem.hh
generallensspatialparams.hh
test_generalproblem2p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/generalproblem)
