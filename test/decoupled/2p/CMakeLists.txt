add_dumux_test(test_impes test_impes test_impes.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_impes-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_impes)

add_dumux_test(test_impesadaptive test_impesadaptive test_impesadaptive.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_2padaptive-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2padaptive-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_impesadaptive)

add_executable("test_impeswithamg" EXCLUDE_FROM_ALL test_impeswithamg.cc)
if(MPI_FOUND)
  add_dumux_test(test_impeswithamg test_impeswithamg test_impeswithamg.cc
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference-parallel.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-test_impeswithamg-00094.vtu
    "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg"
    -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_impeswithamg.input
    -TimeManager.TEnd 7e7)
else()
  add_dumux_test(test_impeswithamg test_impeswithamg test_impeswithamg.cc
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg-00009.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg)
endif(MPI_FOUND)

add_dumux_test(test_transport test_transport test_transport.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_transport-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_transport-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_transport)

add_dumux_test(test_mpfao2p test_mpfa2p test_mpfa2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_mpfao2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input
  -ModelType MPFAO)

add_dumux_test(test_mpfal2p test_mpfa2p test_mpfa2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_mpfal2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input
  -ModelType MPFAL)

add_dumux_test(test_mpfal2padaptive test_mpfa2p test_mpfa2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_mpfal2padaptive-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input
  -ModelType MPFALAdaptive)

add_dumux_test(test_3d2p test_3d2p test_3d2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_3d2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_3d2p-00012.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_3d2p)
