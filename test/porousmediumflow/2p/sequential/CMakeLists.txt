add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_impes
              LABELS porousmediumflow 2p sequential
              SOURCES test_impes.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_impes-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_impes")

dumux_add_test(NAME test_impesadaptive
              LABELS porousmediumflow 2p sequential
              SOURCES test_impesadaptive.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2padaptive-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2padaptive-00006.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_impesadaptive")

dumux_add_test(NAME test_impesadaptiverestart
              LABELS porousmediumflow 2p sequential
              SOURCES test_impesadaptiverestart.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2padaptive-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2padaptive-00006.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_impesadaptiverestart")

# the restart test has to run after the test that produces the restart file
set_tests_properties(test_impesadaptiverestart PROPERTIES DEPENDS test_impesadaptive)

if(MPI_FOUND)
  dumux_add_test(NAME test_impeswithamg
                LABELS porousmediumflow 2p parallel sequential
                SOURCES test_impeswithamg.cc
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMD_ARGS --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference-parallel.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-test_impeswithamg-00093.vtu
                         --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_impeswithamg.input -TimeManager.TEnd 7e7")
else()
  dumux_add_test(NAME test_impeswithamg
                LABELS porousmediumflow 2p sequential
                SOURCES test_impeswithamg.cc
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMD_ARGS --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg-00008.vtu
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg")
endif()

dumux_add_test(NAME test_transport
              LABELS porousmediumflow 2p sequential
              SOURCES test_transport.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_transport-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_transport-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_transport")

# mpfa tests
add_executable(test_mpfa2p test_mpfa2p.cc)
dumux_add_test(NAME test_mpfao2p
              TARGET test_mpfa2p
              LABELS porousmediumflow 2p sequential
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_mpfao2p-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_mpfao2p-00006.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input -Problem.Name test_mpfao2p -ModelType MPFAO")

dumux_add_test(NAME test_mpfal2p
              TARGET test_mpfa2p
              LABELS porousmediumflow 2p sequential
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_mpfal2p-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_mpfal2p-00006.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input -Problem.Name test_mpfal2p -ModelType MPFAL")

dumux_add_test(NAME test_mpfal2padaptive
              TARGET test_mpfa2p
              LABELS porousmediumflow 2p sequential
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_mpfal2padaptive-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_mpfal2padaptive-00006.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input -Problem.Name test_mpfal2padaptive -ModelType MPFALAdaptive")

# 3d tests
dumux_add_test(NAME test_3d2pfv
              LABELS porousmediumflow 2p sequential
              SOURCES test_3d2p.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pfv-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pfv-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pfv -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pfv -ModelType FV")

dumux_add_test(NAME test_3d2pfvadaptive
              LABELS porousmediumflow 2p sequential
              SOURCES test_3d2p.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pfvadaptive-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pfvadaptive-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pfvadaptive -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pfvadaptive -ModelType FVAdaptive")

dumux_add_test(NAME test_3d2pmimetic
              LABELS porousmediumflow 2p sequential
              SOURCES test_3d2p.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pmimetic-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmimetic-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmimetic -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pmimetic -ModelType Mimetic")

dumux_add_test(NAME test_3d2pmimeticadaptive
              LABELS porousmediumflow 2p sequential
              SOURCES test_3d2p.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pmimeticadaptive-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmimeticadaptive-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmimeticadaptive -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pmimeticadaptive -ModelType MimeticAdaptive")

dumux_add_test(NAME test_3d2pmpfal
              LABELS porousmediumflow 2p sequential
              SOURCES test_3d2p.cc
              TIMEOUT 1500
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pmpfal-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmpfal-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmpfal -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pmpfal -ModelType MPFAL")

dumux_add_test(NAME test_3d2pmpfaladaptive
              LABELS porousmediumflow 2p sequential
              SOURCES test_3d2p.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3d2pmpfaladaptive-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmpfaladaptive-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3d2pmpfaladaptive -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_3d2p.input -Problem.OutputName test_3d2pmpfaladaptive -ModelType MPFALAdaptive")
