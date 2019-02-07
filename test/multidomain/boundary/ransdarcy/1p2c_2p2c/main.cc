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
 *
 * \brief A test problem for the nonisothermal compositional coupled RANS/Darcy problem (1p2cni/2p2cni)
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/partial.hh>
#include <dumux/common/dumuxmessage.hh>
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
#include "problem_rans.hh"

namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, Properties::TTag::RANSTYPETAG>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyTwoPTwoCNI>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyTwoPTwoCNI>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::RANSTYPETAG, Properties::TTag::RANSTYPETAG, TypeTag>;
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
    using RANSTypeTag = Properties::TTag::RANSTYPETAG;
    using DarcyTypeTag = Properties::TTag::DarcyTwoPTwoCNI;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    using RANSGridManager = Dumux::GridManager<GetPropType<RANSTypeTag, Properties::Grid>>;
    RANSGridManager ransGridManager;
    ransGridManager.init("RANS"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();
    const auto& ransGridView = ransGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using RANSFVGridGeometry = GetPropType<RANSTypeTag, Properties::FVGridGeometry>;
    auto ransFvGridGeometry = std::make_shared<RANSFVGridGeometry>(ransGridView);
    ransFvGridGeometry->update();
    using DarcyFVGridGeometry = GetPropType<DarcyTypeTag, Properties::FVGridGeometry>;
    auto darcyFvGridGeometry = std::make_shared<DarcyFVGridGeometry>(darcyGridView);
    darcyFvGridGeometry->update();

    using Traits = StaggeredMultiDomainTraits<RANSTypeTag, RANSTypeTag, DarcyTypeTag>;

    // the coupling manager
    using CouplingManager = StokesDarcyCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(ransFvGridGeometry, darcyFvGridGeometry);

    // the indices
    constexpr auto ransCellCenterIdx = CouplingManager::stokesCellCenterIdx;
    constexpr auto ransFaceIdx = CouplingManager::stokesFaceIdx;
    constexpr auto darcyIdx = CouplingManager::darcyIdx;

    // initialize the fluidsystem (tabulation)
    GetPropType<RANSTypeTag, Properties::FluidSystem>::init();

    // the problem (initial and boundary conditions)
    using RANSProblem = GetPropType<RANSTypeTag, Properties::Problem>;
    auto ransProblem = std::make_shared<RANSProblem>(ransFvGridGeometry, couplingManager);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyFvGridGeometry, couplingManager);

    // get some time loop parameters
    using Scalar = GetPropType<RANSTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // set timeloop for the subproblems, needed for boundary value variations
    ransProblem->setTimeLoop(timeLoop);
    darcyProblem->setTimeLoop(timeLoop);

    // the solution vector
    Traits::SolutionVector sol;
    sol[ransCellCenterIdx].resize(ransFvGridGeometry->numCellCenterDofs());
    sol[ransFaceIdx].resize(ransFvGridGeometry->numFaceDofs());
    sol[darcyIdx].resize(darcyFvGridGeometry->numDofs());

    // get a solution vector storing references to the two Stokes solution vectors
    auto ransSol = partial(sol, ransCellCenterIdx, ransFaceIdx);

    // apply initial solution for instationary problems
    ransProblem->applyInitialSolution(ransSol);
    ransProblem->updateStaticWallProperties();
    ransProblem->updateDynamicWallProperties(ransSol);
    darcyProblem->applyInitialSolution(sol[darcyIdx]);

    auto solOld = sol;

    couplingManager->init(ransProblem, darcyProblem, sol);

    // the grid variables
    using RANSGridVariables = GetPropType<RANSTypeTag, Properties::GridVariables>;
    auto ransGridVariables = std::make_shared<RANSGridVariables>(ransProblem, ransFvGridGeometry);
    ransGridVariables->init(ransSol);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyFvGridGeometry);
    darcyGridVariables->init(sol[darcyIdx]);
    darcyProblem->init(sol[darcyIdx], *darcyGridVariables, ransProblem->ransName());

    // intialize the vtk output module
    StaggeredVtkOutputModule<RANSGridVariables, decltype(ransSol)> ransVtkWriter(*ransGridVariables, ransSol, ransProblem->name());
    GetPropType<RANSTypeTag, Properties::IOFields>::initOutputModule(ransVtkWriter);
    ransVtkWriter.write(0.0);

    VtkOutputModule<DarcyGridVariables, GetPropType<DarcyTypeTag, Properties::SolutionVector>> darcyVtkWriter(*darcyGridVariables, sol[darcyIdx], darcyProblem->name());
    using DarcyVelocityOutput = GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
    darcyVtkWriter.addVelocityOutput(std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(darcyVtkWriter);
    darcyVtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(ransProblem, ransProblem, darcyProblem),
                                                 std::make_tuple(ransFvGridGeometry->cellCenterFVGridGeometryPtr(),
                                                                 ransFvGridGeometry->faceFVGridGeometryPtr(),
                                                                 darcyFvGridGeometry),
                                                 std::make_tuple(ransGridVariables->cellCenterGridVariablesPtr(),
                                                                 ransGridVariables->faceGridVariablesPtr(),
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

        // post time step treatment of Darcy problem
        darcyProblem->postTimeStep(sol[darcyIdx], *darcyGridVariables);

        // advance grid variables to the next time step
        ransGridVariables->advanceTimeStep();
        darcyGridVariables->advanceTimeStep();

        ransSol[ransCellCenterIdx] = sol[ransCellCenterIdx];
        ransSol[ransFaceIdx] = sol[ransFaceIdx];

        // update the auxiliary free flow solution vector
        ransProblem->updateDynamicWallProperties(ransSol);
        assembler->updateGridVariables(sol);

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();


        // write vtk output
        ransVtkWriter.write(timeLoop->time());
        darcyVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(ransGridView.comm());
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
