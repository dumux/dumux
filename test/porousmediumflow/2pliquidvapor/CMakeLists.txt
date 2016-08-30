#add links to input files
add_input_file_links()

add_dumux_test(validation_ochs validation_ochs validationochs.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p")

#add executables concept_preheat validation_ochs field_preheat
# add_executable_all(concept_superheated conceptsuperheated.cc)
# add_executable_all(validation_ochs validationochs.cc)
