// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief A test for the shallow water model (lake at rest).
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include "properties.hh"

//! Deviation from the lake-at-rest solution, measured in the maximum and in the L1 norm
struct WellBalancedError
{
    double freeSurfaceMax = 0.0;
    double freeSurfaceL1 = 0.0;
    double dischargeMax = 0.0;
    double dischargeL1 = 0.0;
};

/*!
 * \brief Compute the deviation from the lake-at-rest solution
 *
 * The bed does not move, so the deviation of the water depth from its exact value is also
 * the deviation of the free surface elevation. The exact specific discharge is zero.
 */
template<class Problem, class SolutionVector>
WellBalancedError computeWellBalancedError(const Problem& problem, const SolutionVector& sol)
{
    const auto& gridGeometry = problem.gridGeometry();
    WellBalancedError error;

    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {
        fvGeometry.bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& priVars = sol[scv.dofIndex()];
            const auto waterDepth = priVars[0];

            using std::abs; using std::max; using std::sqrt;
            const auto freeSurfaceError = abs(waterDepth - problem.exactWaterDepthAtPos(scv.center()));
            const auto dischargeError = waterDepth*sqrt(priVars[1]*priVars[1] + priVars[2]*priVars[2]);

            error.freeSurfaceMax = max(error.freeSurfaceMax, freeSurfaceError);
            error.dischargeMax = max(error.dischargeMax, dischargeError);
            error.freeSurfaceL1 += freeSurfaceError*scv.volume();
            error.dischargeL1 += dischargeError*scv.volume();
        }
    }

    const auto& comm = gridGeometry.gridView().comm();
    error.freeSurfaceMax = comm.max(error.freeSurfaceMax);
    error.dischargeMax = comm.max(error.dischargeMax);
    error.freeSurfaceL1 = comm.sum(error.freeSurfaceL1);
    error.dischargeL1 = comm.sum(error.dischargeL1);

    return error;
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::LakeAtRest;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(dt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance the time loop to the next step
        timeLoop->advanceTimeStep();
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());
    vtkWriter.write(timeLoop->time());

    ////////////////////////////////////////////////////////////
    // check that the lake at rest has been preserved
    ////////////////////////////////////////////////////////////

    const auto error = computeWellBalancedError(*problem, x);
    if (leafGridView.comm().rank() == 0)
        std::cout << "Deviation from the lake at rest after " << timeLoop->time() << " seconds:\n"
                  << std::scientific << std::setprecision(6)
                  << "  free surface elevation: " << error.freeSurfaceMax << " m (max), "
                  << error.freeSurfaceL1 << " m^3 (L1)\n"
                  << "  specific discharge:     " << error.dischargeMax << " m^2/s (max), "
                  << error.dischargeL1 << " m^4/s (L1)" << std::endl;

    // a scheme that is not well-balanced produces an error of the order of the bed variation
    const auto tolerance = getParam<Scalar>("Problem.ErrorTolerance", 1e-10);
    if (error.freeSurfaceMax > tolerance || error.dischargeMax > tolerance)
        DUNE_THROW(Dune::Exception, "Lake at rest is not preserved: free surface elevation is off by "
                                    << error.freeSurfaceMax << " m and specific discharge by "
                                    << error.dischargeMax << " m^2/s, tolerance is " << tolerance);

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
