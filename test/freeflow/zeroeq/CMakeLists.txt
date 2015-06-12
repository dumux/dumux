# test: test_zeroeq
add_dumux_test(test_zeroeq test_zeroeq test_zeroeq.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/zeroeq-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/zeroeq-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_zeroeq.input
  )

# test: test_zeroeq_channel_channel
add_dumux_test(test_zeroeq_channel test_zeroeq_channel test_zeroeq_channel.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/zeroeqchannel-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/zeroeqchannel-00022.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq_channel
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_zeroeq_channel.input
  )
