# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_1p_layered_box_2d
               SOURCES main.cc
               LABELS porousmediumflow 1p
               COMPILE_DEFINITIONS TYPETAG=OnePLayeredBox
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
               COMMAND ./test_1p_layered_box_2d
               CMD_ARGS params.input)

dumux_add_test(NAME test_1p_layered_box_3d
               SOURCES main.cc
               LABELS porousmediumflow 1p
               COMPILE_DEFINITIONS TYPETAG=OnePLayeredBox
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<3>
               COMMAND ./test_1p_layered_box_3d
               CMD_ARGS params.input -Grid.UpperRight "1 1 1" -Grid.Cells "8 8 8")

dumux_add_test(NAME test_1p_layered_pq1bubble
               SOURCES main.cc
               LABELS porousmediumflow 1p
               COMPILE_DEFINITIONS TYPETAG=OnePLayeredPQ1Bubble
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
               COMMAND ./test_1p_layered_pq1bubble
               CMD_ARGS params.input)

dumux_add_test(NAME test_1p_layered_pq1bubble_hybrid
               SOURCES main.cc
               LABELS porousmediumflow 1p
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS TYPETAG=OnePLayeredPQ1BubbleHybrid
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMMAND ./test_1p_layered_pq1bubble_hybrid
               CMD_ARGS params.input)
