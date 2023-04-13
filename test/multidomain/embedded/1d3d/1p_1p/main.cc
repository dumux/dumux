// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Test for the 1d-3d embedded mixed-dimension model coupling two
 *        one-phase porous medium flow problems.
 */

#include <config.h>

#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/geometry/diameter.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "properties.hh"
#include "l2norm.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using BulkTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimTypeTag = Properties::TTag::LOWDIMTYPETAG;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using BulkGridManager = Dumux::GridManager<GetPropType<BulkTypeTag, Properties::Grid>>;
    BulkGridManager bulkGridManager;
    bulkGridManager.init("Tissue"); // pass parameter group

    using LowDimGridManager = Dumux::GridManager<GetPropType<LowDimTypeTag, Properties::Grid>>;
    LowDimGridManager lowDimGridManager;
    lowDimGridManager.init("Vessel"); // pass parameter group

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& bulkGridView = bulkGridManager.grid().leafGridView();
    const auto& lowDimGridView = lowDimGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using BulkFVGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    using LowDimFVGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();

    // the coupling manager
    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(bulkFvGridGeometry, lowDimFvGridGeometry);

    // the problem (initial and boundary conditions)
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, couplingManager);
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkFvGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(sol[bulkIdx]);
    lowDimProblem->applyInitialSolution(sol[lowDimIdx]);
    auto oldSol = sol;

    couplingManager->init(bulkProblem, lowDimProblem, sol);
    bulkProblem->computePointSourceMap();
    lowDimProblem->computePointSourceMap();

    // the grid variables
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    bulkGridVariables->init(sol[bulkIdx]);
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    lowDimGridVariables->init(sol[lowDimIdx]);

    const bool writeVtk = getParam<bool>("Vtk.EnableVtkOutput", true);

    // initialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(bulkVtkWriter);
    bulkProblem->addVtkOutputFields(bulkVtkWriter);

    using LowDimSolutionVector = std::decay_t<decltype(sol[lowDimIdx])>;
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
    GetPropType<LowDimTypeTag, Properties::IOFields>::initOutputModule(lowDimVtkWriter);
    lowDimProblem->addVtkOutputFields(lowDimVtkWriter);

    if (writeVtk)
    {
        bulkVtkWriter.write(0.0);
        lowDimVtkWriter.write(0.0);
    }

    // compute hmax of both domains
    double hMaxBulk = 0.0;
    for (const auto& element : elements(bulkGridView))
        hMaxBulk = std::max(hMaxBulk, diameter(element.geometry()));

    double hMaxLowDim = 0.0;
    for (const auto& element : elements(lowDimGridView))
        hMaxLowDim = std::max(hMaxLowDim, diameter(element.geometry()));

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager);

    // the linear solver
    using LinearSolver = BlockDiagILU0BiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    Dune::Timer assembleTimer(false), solveTimer(false), updateTimer(false);
    std::cout << "\nAssembling linear system... " << std::flush;

    // assemble stiffness matrix
    assembleTimer.start();
    couplingManager->updateSolution(sol);
    assembler->assembleJacobianAndResidual(sol);
    assembleTimer.stop();

    std::cout << "done.\n";
    std::cout << "Solving linear system ("
              << linearSolver->name() << ") ... " << std::flush;

    // solve linear system
    solveTimer.start();
    auto deltaSol = sol;
    const bool converged = linearSolver->solve(assembler->jacobian(), deltaSol, assembler->residual());
    if (!converged) DUNE_THROW(Dune::MathError, "Linear solver did not converge!");
    solveTimer.stop();

    // update variables
    updateTimer.start();
    sol -= deltaSol;
    couplingManager->updateSolution(sol);
    assembler->updateGridVariables(sol);
    updateTimer.stop();

    std::cout << "done.\n";
    const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
    std::cout << "Assemble/solve/update time: "
              <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
              <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
              <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
              << "\n";

    // output the source terms
    bulkProblem->computeSourceIntegral(sol[bulkIdx], *bulkGridVariables);
    auto sourceNorm = lowDimProblem->computeSourceIntegral(sol[lowDimIdx], *lowDimGridVariables);

    // write vtk output
    if (writeVtk)
    {
        bulkVtkWriter.write(1.0);
        lowDimVtkWriter.write(1.0);
    }

    // write the L2-norm to file
    static const int order = getParam<int>("Problem.NormIntegrationOrder");
    auto norm1D = L2Norm<double>::computeErrorNorm(*lowDimProblem, sol[lowDimIdx], order);
    norm1D /= L2Norm<double>::computeNormalization(*lowDimProblem, order);
    auto norm3D = L2Norm<double>::computeErrorNorm(*bulkProblem, sol[bulkIdx], order);
    norm3D /= L2Norm<double>::computeNormalization(*bulkProblem, order);

    // output result to terminal
    std::cout << "-----------------------------------------------------\n";
    std::cout << " L2_p1d: " << norm1D << " hmax_1d: " << hMaxLowDim << '\n'
              << " L2_p3d: " << norm3D << " hmax_3d: " << hMaxBulk << '\n'
              << " L2_q  : " << sourceNorm << " hmax_3d: " << hMaxBulk << '\n';
    std::cout << "-----------------------------------------------------" << std::endl;

    // ... and file.
    {
        std::ofstream outFile(getParam<std::string>("Problem.OutputFilename"));
        outFile << hMaxBulk << " " << hMaxLowDim << " " << norm3D << " " << norm1D << " " << sourceNorm << '\n';
    }

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
