add_input_file_links()

add_dumux_test(test_zeroeq test_zeroeq test_zeroeq.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq"
                  --files ${CMAKE_SOURCE_DIR}/test/references/zeroeq-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/zeroeq-00005.vtu)

add_dumux_test(test_zeroeq_channel test_zeroeq_channel test_zeroeq_channel.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq_channel"
                  --files ${CMAKE_SOURCE_DIR}/test/references/zeroeqchannel-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/zeroeqchannel-00021.vtu)

#install sources
install(FILES
test_zeroeq.cc
test_zeroeq_channel.cc
zeroeqchanneltestproblem.hh
zeroeqtestproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/zeroeq)
