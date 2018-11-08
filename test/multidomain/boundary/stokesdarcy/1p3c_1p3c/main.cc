// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief A test problem for the coupled Stokes/Darcy problem (1p)
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/geometry/diameter.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/multidomain/staggeredtraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingmanager.hh>

#include "problem_darcy.hh"
#include "problem_stokes.hh"

namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::StokesOnePThreeCTypeTag>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyOnePThreeCTypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOnePThreeCTypeTag>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::StokesOnePThreeCTypeTag, Properties::TTag::StokesOnePThreeCTypeTag, TypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv) try
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
    using StokesTypeTag = Properties::TTag::StokesOnePThreeCTypeTag;
    using DarcyTypeTag = Properties::TTag::DarcyOnePThreeCTypeTag;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    using StokesGridManager = Dumux::GridManager<GetPropType<StokesTypeTag, Properties::Grid>>;
    StokesGridManager stokesGridManager;
    stokesGridManager.init("Stokes"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();
    const auto& stokesGridView = stokesGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using StokesFVGridGeometry = GetPropType<StokesTypeTag, Properties::FVGridGeometry>;
    auto stokesFvGridGeometry = std::make_shared<StokesFVGridGeometry>(stokesGridView);
    stokesFvGridGeometry->update();
    using DarcyFVGridGeometry = GetPropType<DarcyTypeTag, Properties::FVGridGeometry>;
    auto darcyFvGridGeometry = std::make_shared<DarcyFVGridGeometry>(darcyGridView);
    darcyFvGridGeometry->update();

    using Traits = StaggeredMultiDomainTraits<StokesTypeTag, StokesTypeTag, DarcyTypeTag>;

    // the coupling manager
    using CouplingManager = StokesDarcyCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(stokesFvGridGeometry, darcyFvGridGeometry);

    // the indices
    constexpr auto stokesCellCenterIdx = CouplingManager::stokesCellCenterIdx;
    constexpr auto stokesFaceIdx = CouplingManager::stokesFaceIdx;
    constexpr auto darcyIdx = CouplingManager::darcyIdx;

    // the problems (initial and boundary conditions)
    using StokesProblem = GetPropType<StokesTypeTag, Properties::Problem>;
    auto stokesProblem = std::make_shared<StokesProblem>(stokesFvGridGeometry, couplingManager);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyFvGridGeometry, couplingManager);

    // get some time loop parameters
    using Scalar = GetPropType<StokesTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the solution vector
    Traits::SolutionVector sol;
    sol[stokesCellCenterIdx].resize(stokesFvGridGeometry->numCellCenterDofs());
    sol[stokesFaceIdx].resize(stokesFvGridGeometry->numFaceDofs());
    sol[darcyIdx].resize(darcyFvGridGeometry->numDofs());

    const auto& cellCenterSol = sol[stokesCellCenterIdx];
    const auto& faceSol = sol[stokesFaceIdx];

    // apply initial solution for instationary problems
    GetPropType<StokesTypeTag, Properties::SolutionVector> stokesSol;
    std::get<0>(stokesSol) = cellCenterSol;
    std::get<1>(stokesSol) = faceSol;
    stokesProblem->applyInitialSolution(stokesSol);
    sol[stokesCellCenterIdx] = stokesSol[stokesCellCenterIdx];
    sol[stokesFaceIdx] = stokesSol[stokesFaceIdx];

    darcyProblem->applyInitialSolution(sol[darcyIdx]);
    auto solDarcyOld = sol[darcyIdx];

    auto solOld = sol;

    couplingManager->init(stokesProblem, darcyProblem, sol);

    // the grid variables
    using StokesGridVariables = GetPropType<StokesTypeTag, Properties::GridVariables>;
    auto stokesGridVariables = std::make_shared<StokesGridVariables>(stokesProblem, stokesFvGridGeometry);
    stokesGridVariables->init(stokesSol);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyFvGridGeometry);
    darcyGridVariables->init(sol[darcyIdx]);

    // intialize the vtk output module
    const auto stokesName = getParam<std::string>("Problem.Name") + "_" + stokesProblem->name();
    const auto darcyName = getParam<std::string>("Problem.Name") + "_" + darcyProblem->name();

    StaggeredVtkOutputModule<StokesGridVariables, GetPropType<StokesTypeTag, Properties::SolutionVector>> stokesVtkWriter(*stokesGridVariables, stokesSol, stokesName);
    GetPropType<StokesTypeTag, Properties::IOFields>::initOutputModule(stokesVtkWriter);
    stokesVtkWriter.write(0.0);

    VtkOutputModule<DarcyGridVariables, GetPropType<DarcyTypeTag, Properties::SolutionVector>> darcyVtkWriter(*darcyGridVariables, sol[darcyIdx], darcyName);
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(darcyVtkWriter);
    darcyVtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(stokesProblem, stokesProblem, darcyProblem),
                                                 std::make_tuple(stokesFvGridGeometry->cellCenterFVGridGeometryPtr(),
                                                                 stokesFvGridGeometry->faceFVGridGeometryPtr(),
                                                                 darcyFvGridGeometry),
                                                 std::make_tuple(stokesGridVariables->cellCenterGridVariablesPtr(),
                                                                 stokesGridVariables->faceGridVariablesPtr(),
                                                                 darcyGridVariables),
                                                 couplingManager,
                                                 timeLoop);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(solOld);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        solOld = sol;
        stokesGridVariables->advanceTimeStep();
        darcyGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        stokesVtkWriter.write(timeLoop->time());
        darcyVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(stokesGridView.comm());
    timeLoop->finalize(darcyGridView.comm());

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
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
