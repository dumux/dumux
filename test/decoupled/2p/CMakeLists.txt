# add build targets
add_executable("test_impes" EXCLUDE_FROM_ALL test_impes.cc)
add_test(test_impes
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_impes-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_impes
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_impes.input)

add_executable("test_impesadaptive" EXCLUDE_FROM_ALL test_impesadaptive.cc)
add_dune_alugrid_flags(test_impesadaptive)
add_test(test_impesadaptive
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_2padaptive-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2padaptive-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_impesadaptive
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_impesadaptive.input)

add_executable("test_impeswithamg" EXCLUDE_FROM_ALL test_impeswithamg.cc)
if(MPI_FOUND)
  add_dune_superlu_flags(test_impeswithamg)
  add_dune_mpi_flags(test_impeswithamg)
  add_test(test_impeswithamg
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference-parallel.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-test_impeswithamg-00013.vtu
    "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg"
    -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_impesadaptive.input
    -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_impeswithamg.dgf)
else(MPI_FOUND)
  add_test(test_impeswithamg
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/test_impes-reference.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg-00013.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/test_impeswithamg
    -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_impesadaptive.input
    -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_impeswithamg.dgf)
endif(MPI_FOUND)

add_executable("test_transport" EXCLUDE_FROM_ALL test_transport.cc)
add_test(test_transport
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_transport-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_transport-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_transport
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_transport.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_transport.dgf)

add_executable("test_mpfa2p" EXCLUDE_FROM_ALL test_mpfa2p.cc)
add_dune_alugrid_flags(test_mpfa2p)
add_dune_superlu_flags(test_mpfa2p)
add_test(test_mpfao2p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_mpfao2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input
  -ModelType MPFAO)
add_test(test_mpfal2p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_mpfal2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input
  -ModelType MPFAL)
add_test(test_mpfal2padaptive
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_mpfal2padaptive-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_mpfa2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_mpfa2p.input
  -ModelType MPFALAdaptive)

set(NORMALTESTS
  test_impes
  test_impesadaptive
  test_impeswithamg
  test_transport
  test_mpfa2p)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY TEST_INCLUDE_FILE ${CMAKE_CURRENT_BINARY_DIR}/BuildTests.cmake)

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
