// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief Test for the equal dimension boundary coupling model.
 */

#include <config.h>

#include <ctime>
#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/boundary/darcydarcy/couplingmanager.hh>

#include "problem.hh"

#ifndef DOMAINSPLIT
#define DOMAINSPLIT 0
#endif

namespace Dumux {
namespace Properties {

// Create new type tags
namespace TTag {
struct OnePSub { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
// differentiate between the two subproblems
struct OnePSub0 { using InheritsFrom = std::tuple<OnePSub>; };
struct OnePSub1 { using InheritsFrom = std::tuple<OnePSub>; };
} // end namespace TTag

// the coupling manager
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePSub>
{ using type = DarcyDarcyBoundaryCouplingManager<MultiDomainTraits<Properties::TTag::OnePSub0, Properties::TTag::OnePSub1>>; };

// Set the grid type
#if DOMAINSPLIT==1
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePSub>
{
    using HostGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
    using type = Dune::SubGrid<HostGrid::dimension, HostGrid>;
};
#elif DOMAINSPLIT==0
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePSub>
{
    using HostGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
    using type = HostGrid;
};
#endif

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePSub>
{
    using type = OnePTestSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                       GetPropType<TypeTag, Properties::Scalar>>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePSub>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// differentiate between the two subproblems
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePSub0> { using type = OnePTestProblem<TypeTag, 0>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePSub1> { using type = OnePTestProblem<TypeTag, 1>; };

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using TypeTag = Properties::TTag::OnePSub;
    using SubTypeTag0 = Properties::TTag::OnePSub0;
    using SubTypeTag1 = Properties::TTag::OnePSub1;

    // create the full grid that we are gonna split for the output
    using HostGrid = typename GetProp<TypeTag, Properties::Grid>::HostGrid;
    GridManager<HostGrid> gridManager;
    gridManager.init();

    // get the lens
    using GlobalPosition = typename HostGrid::template Codim<0>::Geometry::GlobalCoordinate;
    const auto lensLowerLeft = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
    const auto lensUpperRight = getParam<GlobalPosition>("SpatialParams.LensUpperRight");

    // grid managers for the subdomains
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager0, gridManager1;

///////////////////////////////////////////////////////////////////////////////////////////
// split the domains by creating two separate grids for lens and the rest (using sub-grid)
///////////////////////////////////////////////////////////////////////////////////////////
#if DOMAINSPLIT==1

    // create subgrids coinciding with the lens
    auto elementSelector0 = [&lensLowerLeft, &lensUpperRight](const auto& element)
    { return LensSpatialParams::pointInLens(element.geometry().center(), lensLowerLeft, lensUpperRight); };
    auto elementSelector1 = [&lensLowerLeft, &lensUpperRight](const auto& element)
    { return !LensSpatialParams::pointInLens(element.geometry().center(), lensLowerLeft, lensUpperRight); };

    gridManager0.init(gridManager.grid(), elementSelector0, "1");
    gridManager1.init(gridManager.grid(), elementSelector1, "2");

//////////////////////////////////////////////////////////////////////////////////
// split the domains by creating two separate grids for the lower and upper half
//////////////////////////////////////////////////////////////////////////////////
#elif DOMAINSPLIT==0

    // create an upper half and lower half grid
    gridManager0.init("1");
    gridManager1.init("2");

#endif

    // we compute on the leaf grid views
    const auto& gridView0 = gridManager0.grid().leafGridView();
    const auto& gridView1 = gridManager1.grid().leafGridView();

    ////////////////////////////////////////////////
    // run the multidomain simulation on two grids
    ////////////////////////////////////////////////

    // create the finite volume grid geometries
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto fvGridGeometry0 = std::make_shared<GridGeometry>(gridView0);
    auto fvGridGeometry1 = std::make_shared<GridGeometry>(gridView1);
    fvGridGeometry0->update();
    fvGridGeometry1->update();

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<SubTypeTag0, SubTypeTag1>;
    constexpr auto domain0Idx = Traits::template SubDomain<0>::Index();
    constexpr auto domain1Idx = Traits::template SubDomain<1>::Index();

    // the coupling manager
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (initial and boundary conditions)
    using Problem0 = GetPropType<SubTypeTag0, Properties::Problem>;
    auto problem0 = std::make_shared<Problem0>(fvGridGeometry0, couplingManager, "1");
    problem0->spatialParams().setLens(lensLowerLeft, lensUpperRight);
    using Problem1 = GetPropType<SubTypeTag1, Properties::Problem>;
    auto problem1 = std::make_shared<Problem1>(fvGridGeometry1, couplingManager, "2");
    problem1->spatialParams().setLens(lensLowerLeft, lensUpperRight);

    // the solution vector
    Traits::SolutionVector sol;
    sol[domain0Idx].resize(fvGridGeometry0->numDofs());
    sol[domain1Idx].resize(fvGridGeometry1->numDofs());
    problem0->applyInitialSolution(sol[domain0Idx]);
    problem1->applyInitialSolution(sol[domain1Idx]);
    auto oldSol = sol;

    // compute the coupling map
    couplingManager->init(problem0, problem1, sol);

    // the grid variables
    using GridVariables0 = GetPropType<SubTypeTag0, Properties::GridVariables>;
    auto gridVariables0 = std::make_shared<GridVariables0>(problem0, fvGridGeometry0);
    using GridVariables1 = GetPropType<SubTypeTag1, Properties::GridVariables>;
    auto gridVariables1 = std::make_shared<GridVariables1>(problem1, fvGridGeometry1);
    gridVariables0->init(sol[domain0Idx]);
    gridVariables1->init(sol[domain1Idx]);

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // intialize the vtk output module
    using SolutionVector0 = std::decay_t<decltype(sol[domain0Idx])>;
    VtkOutputModule<GridVariables0, SolutionVector0> vtkWriter0(*gridVariables0, sol[domain0Idx], problem0->name());
    GetPropType<SubTypeTag0, Properties::IOFields>::initOutputModule(vtkWriter0);
    vtkWriter0.write(0.0);

    using SolutionVector1 = std::decay_t<decltype(sol[domain1Idx])>;
    VtkOutputModule<GridVariables1, SolutionVector1> vtkWriter1(*gridVariables1, sol[domain1Idx], problem1->name());
    GetPropType<SubTypeTag1, Properties::IOFields>::initOutputModule(vtkWriter1);
    vtkWriter1.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(problem0, problem1),
                                                 std::make_tuple(fvGridGeometry0, fvGridGeometry1),
                                                 std::make_tuple(gridVariables0, gridVariables1),
                                                 couplingManager, timeLoop, oldSol);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        oldSol = sol;
        gridVariables0->advanceTimeStep();
        gridVariables1->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        vtkWriter0.write(timeLoop->time());
        vtkWriter1.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }

    timeLoop->finalize(mpiHelper.getCollectiveCommunication());

    //////////////////////////////////////////////////////////////////////////
    // write out a combined vtu for vtu comparison in the testing framework
    //////////////////////////////////////////////////////////////////////////

    const auto& gridView = gridManager.grid().leafGridView();
    CCTpfaFVGridGeometry<typename HostGrid::LeafGridView> gridGeometry(gridView);
    const auto& bBoxTree = gridGeometry.boundingBoxTree();
    // copy data from the subdomains to full domain data vectors
    std::vector<int> processRank(gridView.size(0), 0); // sequential simulation
    std::vector<int> pressure(gridView.size(0));

    for (const auto& element : elements(gridView0))
    {
        const auto eIdx = fvGridGeometry0->elementMapper().index(element);
        const auto eIdxHost = intersectingEntities(element.geometry().center(), bBoxTree)[0];
        pressure[eIdxHost] = sol[domain0Idx][eIdx][0];
    }

    for (const auto& element : elements(gridView1))
    {
        const auto eIdx = fvGridGeometry1->elementMapper().index(element);
        const auto eIdxHost = intersectingEntities(element.geometry().center(), bBoxTree)[0];
        pressure[eIdxHost] = sol[domain1Idx][eIdx][0];
    }

    Dune::VTKWriter<typename HostGrid::LeafGridView> vtkWriter(gridView);
    vtkWriter.addCellData(processRank, "process rank");
    vtkWriter.addCellData(pressure, "pressure");
    const auto filename = getParam<std::string>("Vtk.OutputName") + "_combined";
    vtkWriter.write(filename);

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
