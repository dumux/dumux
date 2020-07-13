dune_add_test(NAME test_throattransmissibility
              SOURCES test_throattransmissibility.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_throattransmissibility)

dune_add_test(NAME test_material_2p_porenetwork
              SOURCES test_material_2p_porenetwork.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_material_2p_porenetwork)
