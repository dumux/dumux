# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_md_embedded_1d3d_parallel_line_tpfatpfa
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Line
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_line_tpfatpfa_bisection
               TARGET test_md_embedded_1d3d_parallel_line_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_line_boxbox
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Line
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_surface_tpfatpfa
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_surface_tpfatpfa_bisection
               TARGET test_md_embedded_1d3d_parallel_surface_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_surface_tpfatpfa_pointvalue
               TARGET test_md_embedded_1d3d_parallel_surface_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection -MixedDimension.UseCircleAverage false)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_surface_boxbox
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_surface_boxbox_nearfaces
               TARGET test_md_embedded_1d3d_parallel_surface_boxbox
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection -Network.NearFaceVessels true)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_average_tpfatpfa
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_average_tpfatpfa_bisection
               TARGET test_md_embedded_1d3d_parallel_average_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_average_boxbox
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Average
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_average_boxbox_nearfaces
               TARGET test_md_embedded_1d3d_parallel_average_boxbox
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection -Network.NearFaceVessels true)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_kernel_tpfatpfa
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_kernel_tpfatpfa_bisection
               TARGET test_md_embedded_1d3d_parallel_kernel_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_kernel_tpfatpfa_order3
               TARGET test_md_embedded_1d3d_parallel_kernel_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection -MixedDimension.KernelIntegrationOrder 3)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_kernel_boxbox
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_line_gradedtpfa
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueGradedCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Line
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_surface_gradedtpfa
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueGradedCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_average_gradedbox
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueGradedBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Average
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_kernel_gradedtpfa
               SOURCES main.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueGradedCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_line_tpfatpfa_emptyprocesses
               TARGET test_md_embedded_1d3d_parallel_line_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.UpperRight "3 1 1" -Tissue.Grid.Cells "30 10 10")

dumux_add_test(NAME test_md_embedded_1d3d_parallel_surface_tpfatpfa_emptyprocesses
               TARGET test_md_embedded_1d3d_parallel_surface_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.UpperRight "3 1 1" -Tissue.Grid.Cells "30 10 10")

dumux_add_test(NAME test_md_embedded_1d3d_parallel_average_tpfatpfa_emptyprocesses
               TARGET test_md_embedded_1d3d_parallel_average_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.UpperRight "3 1 1" -Tissue.Grid.Cells "30 10 10")

dumux_add_test(NAME test_md_embedded_1d3d_parallel_kernel_tpfatpfa_emptyprocesses
               TARGET test_md_embedded_1d3d_parallel_kernel_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.UpperRight "3 1 1" -Tissue.Grid.Cells "30 10 10")

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_line_tpfatpfa
               SOURCES main_jacobian.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Line
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_line_tpfatpfa_bisection
               TARGET test_md_embedded_1d3d_parallel_jacobian_line_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_line_boxbox
               SOURCES main_jacobian.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Line
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_surface_tpfatpfa
               SOURCES main_jacobian.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_surface_tpfatpfa_bisection
               TARGET test_md_embedded_1d3d_parallel_jacobian_surface_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_surface_boxbox
               SOURCES main_jacobian.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_average_tpfatpfa
               SOURCES main_jacobian.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_average_tpfatpfa_bisection
               TARGET test_md_embedded_1d3d_parallel_jacobian_average_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_average_boxbox
               SOURCES main_jacobian.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Average
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_kernel_tpfatpfa
               SOURCES main_jacobian.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_kernel_tpfatpfa_bisection
               TARGET test_md_embedded_1d3d_parallel_jacobian_kernel_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_jacobian_kernel_boxbox
               SOURCES main_jacobian.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_line_tpfatpfa
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Line
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_line_boxbox
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Line
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_tpfatpfa
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_boxbox
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_average_tpfatpfa
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_average_boxbox
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Average
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_kernel_tpfatpfa
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_kernel_boxbox
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_boxtpfa
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_tpfabox
               SOURCES main_solver.cc
               LABELS multidomain multidomain_embedded parallel
               COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkBox COUPLINGMODE=Embedded1d3dCouplingMode::Surface
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Network.Grid.Overlap 0 -Network.Partition Bisection)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_tpfatpfa_gmres
               TARGET test_md_embedded_1d3d_parallel_solver_surface_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Test.Krylov GMRes)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_boxbox_gmres
               TARGET test_md_embedded_1d3d_parallel_solver_surface_boxbox
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection -Test.Krylov GMRes)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_boxtpfa_gmres
               TARGET test_md_embedded_1d3d_parallel_solver_surface_boxtpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Partition Bisection -Test.Krylov GMRes)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_kernel_tpfatpfa_gaussseidel
               TARGET test_md_embedded_1d3d_parallel_solver_kernel_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -LinearSolver.BlockPreconditioner GaussSeidel)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_kernel_boxbox_gaussseidel
               TARGET test_md_embedded_1d3d_parallel_solver_kernel_boxbox
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Grid.Overlap 0 -Network.Partition Bisection -LinearSolver.BlockPreconditioner GaussSeidel)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_boxtpfa_gaussseidel
               TARGET test_md_embedded_1d3d_parallel_solver_surface_boxtpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Tissue.Grid.Overlap 0 -Network.Partition Bisection -LinearSolver.BlockPreconditioner GaussSeidel)

dumux_add_test(NAME test_md_embedded_1d3d_parallel_solver_surface_tpfatpfa_gmres_gaussseidel
               TARGET test_md_embedded_1d3d_parallel_solver_surface_tpfatpfa
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND AND HAVE_UMFPACK )"
               MPI_RANKS 1 2 4
               TIMEOUT 1200
               CMD_ARGS params.input -Test.Krylov GMRes -LinearSolver.BlockPreconditioner GaussSeidel)
