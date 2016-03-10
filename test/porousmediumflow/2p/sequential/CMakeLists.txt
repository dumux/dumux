add_input_file_links()

add_dumux_test(test_impes test_impes test_impes.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_impes-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_impes")

add_dumux_test(test_impesadaptive test_impesadaptive test_impesadaptive.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_2padaptive-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_2padaptive-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_impesadaptive")

add_dumux_test(test_impesadaptiverestart test_impesadaptiverestart test_impesadaptiverestart.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_2padaptive-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_2padaptive-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_impesadaptiverestart")

# the restart test has to run after the test that produces the restart file
set_tests_properties(test_impesadaptiverestart PROPERTIES DEPENDS test_impesadaptive)

if(MPI_FOUND)
  add_dumux_test(test_impeswithamg test_impeswithamg test_impeswithamg.cc
                 python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                   --script fuzzy
                   --files ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference-parallel.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-test_impeswithamg-00094.vtu
                   --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_impeswithamg.input -TimeManager.TEnd 7e7")
else()
  add_dumux_test(test_impeswithamg test_impeswithamg test_impeswithamg.cc
                 python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                   --script fuzzy
                   --files ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg-00009.vtu
                   --command "${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg")
endif()

add_dumux_test(test_transport test_transport test_transport.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_transport-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_transport-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_transport")

add_dumux_test(test_mpfao2p test_mpfa2p test_mpfa2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_mpfao2p-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input -ModelType MPFAO")

add_dumux_test(test_mpfal2p test_mpfa2p test_mpfa2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_mpfal2p-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input -ModelType MPFAL")

add_dumux_test(test_mpfal2padaptive test_mpfa2p test_mpfa2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_mpfal2padaptive-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input -ModelType MPFALAdaptive")

add_dumux_test(test_3d2pfv test_3d2pfv test_3d2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pfv-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pfv-00012.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pfv -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pfv -ModelType FV")

add_dumux_test(test_3d2pfvadaptive test_3d2pfvadaptive test_3d2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pfvadaptive-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pfvadaptive-00012.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pfvadaptive -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pfvadaptive -ModelType FVAdaptive")

add_dumux_test(test_3d2pmimetic test_3d2pmimetic test_3d2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pmimetic-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmimetic-00012.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmimetic -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pmimetic -ModelType Mimetic")

add_dumux_test(test_3d2pmimeticadaptive test_3d2pmimeticadaptive test_3d2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pmimeticadaptive-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmimeticadaptive-00012.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmimeticadaptive -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pmimeticadaptive -ModelType MimeticAdaptive")

add_dumux_test(test_3d2pmpfal test_3d2pmpfal test_3d2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pmpfal-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmpfal-00012.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmpfal -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pmpfal -ModelType MPFAL")

add_dumux_test(test_3d2pmpfaladaptive test_3d2pmpfaladaptive test_3d2p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pmpfaladaptive-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmpfaladaptive-00012.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmpfaladaptive -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pmpfaladaptive -ModelType MPFALAdaptive")

#install sources
install(FILES
buckleyleverettanalyticsolution.hh
mcwhorteranalyticsolution.hh
test_3d2p.cc
test_3d2pproblem.hh
test_3d2pspatialparams.hh
test_impesadaptive.cc
test_impesadaptiveproblem.hh
test_impesadaptivespatialparams.hh
test_impes.cc
test_impesproblem.hh
test_impesspatialparams.hh
test_impeswithamg.cc
test_mpfa2p.cc
test_mpfa2pproblem.hh
test_mpfa2pspatialparams.hh
test_transport.cc
test_transportproblem.hh
test_transportspatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/2p/sequential)
