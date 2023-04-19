// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for the pore-network model with non-creeping flow
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::PNMOnePNonCreepingProblem;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using GridManager = Dumux::PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    // the spatial parameters
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams);

    // the solution vector
    using GridView = typename GridGeometry::GridView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(leafGridView.size(GridView::dimension));
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // the assembler
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // helper class to calculate the boundary flux
    const auto boundaryFlux = Dumux::PoreNetwork::BoundaryFlux(*gridVariables, assembler->localResidual(), x);

    // the linear solver
    using LinearSolver = ILURestartedGMResIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    // output result to vtk
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    Dumux::PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector>
        vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter);

    vtkWriter.write(0.0);

    if (mpiHelper.rank() == 0)
    {
        std::cout << "cumulative outflux is: "
                  << boundaryFlux.getFlux("max", 0, true) << std::endl;
        DumuxMessage::print();
        Parameters::print();
    }

    return 0;
}
