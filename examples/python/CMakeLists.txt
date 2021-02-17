dune_symlink_to_source_files(FILES explicit_transport_cctpfa.py)

dune_python_add_test(NAME explicit_transport_cctpfa
                     SOURCES explicit_transport_cctpfa.py
                     LABELS transport python example
                     CMAKE_GUARD DUNE_ENABLE_PYTHONBINDINGS
                     SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/example-python-explicit-transport-cctpfa-100.vtu
                                     finitevolume-solution-100.vtu
                             --command "${CMAKE_CURRENT_SOURCE_DIR}/explicit_transport_cctpfa.py")
