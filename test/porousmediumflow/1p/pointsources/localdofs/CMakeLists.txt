# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params_1d.input params_2d.input params_2d_simplex.input params_3d.input)

function(add_pointsources_localdofs_test name typetag order gridtype params guard)
    dumux_add_test(NAME ${name}
                   SOURCES main.cc
                   CMAKE_GUARD "${guard}"
                   COMPILE_DEFINITIONS TYPETAG=${typetag}
                   COMPILE_DEFINITIONS ORDER=${order}
                   COMPILE_DEFINITIONS GRIDTYPE=${gridtype}
                   LABELS porousmediumflow 1p
                   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${name}
                   CMD_ARGS ${params} -Problem.Name ${name})
endfunction()

set(umfpack HAVE_UMFPACK)
set(umfpack_ug "( HAVE_UMFPACK AND dune-uggrid_FOUND )")

add_pointsources_localdofs_test(test_1p_pointsources_localdofs_1d_pq1fe OnePPointSourcePQ1FE 1 Dune::OneDGrid params_1d.input "${umfpack}")
add_pointsources_localdofs_test(test_1p_pointsources_localdofs_1d_pq2fe OnePPointSourcePQ2FE 2 Dune::OneDGrid params_1d.input "${umfpack}")

add_pointsources_localdofs_test(test_1p_pointsources_localdofs_2d_box OnePPointSourceBox 1 Dune::YaspGrid<2> params_2d.input "${umfpack}")
add_pointsources_localdofs_test(test_1p_pointsources_localdofs_2d_pq2hybrid OnePPointSourcePQ2Hybrid 2 Dune::YaspGrid<2> params_2d.input "${umfpack}")
add_pointsources_localdofs_test(test_1p_pointsources_localdofs_2d_pq1fe OnePPointSourcePQ1FE 1 Dune::YaspGrid<2> params_2d.input "${umfpack}")
add_pointsources_localdofs_test(test_1p_pointsources_localdofs_2d_pq2fe OnePPointSourcePQ2FE 2 Dune::YaspGrid<2> params_2d.input "${umfpack}")

add_pointsources_localdofs_test(test_1p_pointsources_localdofs_2d_simplex_box OnePPointSourceBox 1 Dune::UGGrid<2> params_2d_simplex.input "${umfpack_ug}")
add_pointsources_localdofs_test(test_1p_pointsources_localdofs_2d_simplex_pq2hybrid OnePPointSourcePQ2Hybrid 2 Dune::UGGrid<2> params_2d_simplex.input "${umfpack_ug}")

add_pointsources_localdofs_test(test_1p_pointsources_localdofs_3d_box OnePPointSourceBox 1 Dune::YaspGrid<3> params_3d.input "${umfpack}")
add_pointsources_localdofs_test(test_1p_pointsources_localdofs_3d_pq2hybrid OnePPointSourcePQ2Hybrid 2 Dune::YaspGrid<3> params_3d.input "${umfpack}")
add_pointsources_localdofs_test(test_1p_pointsources_localdofs_3d_pq2fe OnePPointSourcePQ2FE 2 Dune::YaspGrid<3> params_3d.input "${umfpack}")
