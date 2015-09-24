add_dumux_test(test_box2pmincvol test_box2pminc test_box2pminc.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
               --script fuzzy
               --files ${CMAKE_SOURCE_DIR}/test/references/box2pmincvol-reference.vtu
               ${CMAKE_CURRENT_BINARY_DIR}/box2pmincvol-00027.vtu
               --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pminc -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2pmincvol.input"
               )
  
add_dumux_test(test_box2pmincdist test_box2pminc test_box2pminc.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
               --script fuzzy
               --files ${CMAKE_SOURCE_DIR}/test/references/box2pmincdist-reference.vtu
               ${CMAKE_CURRENT_BINARY_DIR}/box2pmincdist-00029.vtu
               --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pminc -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2pmincdist.input"
               )

# install sources
install(FILES
  2pminctestproblem.hh
  2pminctestspatialparams.hh
  test_box2pminc.cc
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pminc)
