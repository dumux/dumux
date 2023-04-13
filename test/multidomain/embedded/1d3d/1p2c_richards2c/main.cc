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
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/glue.hh>

#include "properties.hh"

namespace Dumux {

//! Helper function for mass balance evaluations
template<class Problem, class SolutionVector, class GridVariables>
double computeSourceIntegral(const Problem& problem, const SolutionVector& sol, const GridVariables& gridVars)
{
    const auto& gg = problem.gridGeometry();
    typename SolutionVector::block_type source(0.0);
    auto fvGeometry = localView(gg);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    for (const auto& element : elements(gg.gridView()))
    {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (auto&& scv : scvs(fvGeometry))
        {
            auto pointSources = problem.scvPointSources(element, fvGeometry, elemVolVars, scv);
            // conversion to kg/s
            const auto& volVars = elemVolVars[scv];
            pointSources *= scv.volume()*volVars.extrusionFactor()
                            * volVars.density(Problem::Indices::liquidPhaseIdx) / volVars.molarDensity(Problem::Indices::liquidPhaseIdx);

            source += pointSources;
        }
    }

    std::cout << "Global integrated source (" << problem.name() << "): " << source[Problem::Indices::conti0EqIdx] << " (kg/s) / "
              <<                           source[Problem::Indices::conti0EqIdx]*3600*24*1000 << " (g/day)" << '\n';

    return source[Problem::Indices::conti0EqIdx];
}

//! Helper function for mass balance evaluations
template<class Problem, class SolutionVector, class GridVariables>
double computeGlobalMass(const Problem& problem, const SolutionVector& sol, const GridVariables& gridVars)
{
    static constexpr int liquidPhaseIdx = Problem::Indices::liquidPhaseIdx;
    static constexpr int transportCompIdx = Problem::Indices::transportCompIdx;
    double mass = 0.0;

    const auto& gg = problem.gridGeometry();
    for (const auto& element : elements(gg.gridView()))
    {
        const auto fvGeometry = localView(gg).bindElement(element);
        const auto elemVolVars = localView(gridVars.curGridVolVars()).bindElement(element, fvGeometry, sol);

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            mass += volVars.massFraction(liquidPhaseIdx, transportCompIdx)*volVars.density(liquidPhaseIdx)
                     *scv.volume() * volVars.porosity() * volVars.saturation(liquidPhaseIdx) * volVars.extrusionFactor();
        }
    }

    return mass;
}

//! Helper function for mass balance evaluations
template<class Problem, class SolutionVector, class GridVariables>
double computeGlobalBoundaryMass(const Problem& problem, const SolutionVector& sol, const GridVariables& gridVars, double dt)
{
    static constexpr int transportCompIdx = Problem::Indices::transportCompIdx;
    static constexpr int transportEqIdx = Problem::Indices::transportEqIdx;
    double mass = 0.0;

    const auto& gg = problem.gridGeometry();
    for (const auto& element : elements(gg.gridView()))
    {
        const auto fvGeometry = localView(gg).bind(element);
        const auto elemVolVars = localView(gridVars.curGridVolVars()).bind(element, fvGeometry, sol);
        const auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache()).bind(element, fvGeometry, elemVolVars);

        for (auto&& scvf : scvfs(fvGeometry))
            if (scvf.boundary())
                mass += problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf)[transportEqIdx]
                        * scvf.area() * elemVolVars[scvf.insideScvIdx()].extrusionFactor()
                        * Problem::FluidSystem::molarMass(transportCompIdx)
                        * dt;
    }

    return mass;
}

} // end namespace Dumux

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
    using BulkTypeTag = Properties::TTag::Soil;
    using LowDimTypeTag = Properties::TTag::Root;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using BulkGridManager = Dumux::GridManager<GetPropType<BulkTypeTag, Properties::Grid>>;
    BulkGridManager bulkGridManager;
    bulkGridManager.init("Soil"); // pass parameter group

    using LowDimGridManager = Dumux::GridManager<GetPropType<LowDimTypeTag, Properties::Grid>>;
    LowDimGridManager lowDimGridManager;
    lowDimGridManager.init("Root"); // pass parameter group

    // we compute on the leaf grid view
    const auto& bulkGridView = bulkGridManager.grid().leafGridView();
    const auto& lowDimGridView = lowDimGridManager.grid().leafGridView();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

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

    // the bulk problem (initial and boundary conditions)
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, couplingManager);

    // the low dim spatial parameters
    using LowDimSpatialParams = GetPropType<LowDimTypeTag, Properties::SpatialParams>;
    auto lowDimSpatialParams = std::make_shared<LowDimSpatialParams>(lowDimFvGridGeometry, lowDimGridManager.getGridData());

    // the low dim problem (initial and boundary conditions)
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager);

    // locally refine levels deep around the embedded grid
    int levels = getParam<int>("Soil.Grid.LocalRefinement");
    for (int i = 0; i < levels; ++i)
    {
        auto& soilGrid = bulkGridManager.grid();
        const auto glue = makeGlue(*lowDimFvGridGeometry, *bulkFvGridGeometry);

        // refine all 3D cells intersected
        for (const auto& is : intersections(glue))
        {
            for (unsigned int targetIdx = 0; targetIdx < is.numTargetNeighbors(); ++targetIdx)
            {
                const auto cutElement = is.targetEntity(targetIdx);

                // mark the cut element and all it's neighbors
                soilGrid.mark(1, cutElement);
                for (const auto& intersection : intersections(bulkGridView, cutElement))
                    if (intersection.neighbor())
                        soilGrid.mark(1, intersection.outside());
            }

        }

        // refine all 3D cells that are where the contamination is
        const double extend = 0.15*(bulkFvGridGeometry->bBoxMax()[0]-bulkFvGridGeometry->bBoxMin()[0]);
        for (const auto& element : elements(bulkGridView))
        {
            const auto globalPos = element.geometry().center();
            auto contaminationPos = bulkFvGridGeometry->bBoxMax()-bulkFvGridGeometry->bBoxMin();
            contaminationPos[0] *= 0.25;
            contaminationPos[1] *= 0.55;
            contaminationPos[2] *= 0.25;
            contaminationPos += bulkFvGridGeometry->bBoxMin();

            if ((globalPos - contaminationPos).infinity_norm() <  extend + 1e-7)
                soilGrid.mark(1, element);
        }

        soilGrid.preAdapt();
        soilGrid.adapt();
        soilGrid.postAdapt();

        // make sure there is only one level difference
        for (int i = 0; i < levels; ++i)
        {
            for (const auto& element : elements(bulkGridView))
            {
                for (const auto& intersection : intersections(bulkGridView, element))
                {
                    if (intersection.neighbor())
                        if (intersection.outside().level()-1 > element.level())
                            soilGrid.mark(1, element);
                }
            }

            soilGrid.preAdapt();
            soilGrid.adapt();
            soilGrid.postAdapt();
        }

        // update the grid geometry (including the bounding box tree)
        bulkFvGridGeometry->update(bulkGridView);
    }

    // output min max h
    double bulkHMin = 1.0; double bulkHMax = 0.0;
    for (const auto& element : elements(bulkGridView))
    {
        const auto geometry = element.geometry();
        const auto h = (geometry.corner(1)-geometry.corner(0)).two_norm();
        bulkHMin = std::min(bulkHMin, h);
        bulkHMax = std::max(bulkHMax, h);
    }

    double ldHMin = 1.0; double ldHMax = 0.0;
    for (const auto& element : elements(lowDimGridView))
    {
        const auto geometry = element.geometry();
        const auto h = (geometry.corner(1)-geometry.corner(0)).two_norm();
        ldHMin = std::min(ldHMin, h);
        ldHMax = std::max(ldHMax, h);
    }

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

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto episodeLength = getParam<Scalar>("TimeLoop.EpisodeLength");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    const bool outputVtk = getParam<bool>("Problem.EnableVtkOutput", true);

    // initialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(bulkVtkWriter);
    bulkVtkWriter.write(0.0);

    using LowDimSolutionVector = std::decay_t<decltype(sol[lowDimIdx])>;
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
    GetPropType<LowDimTypeTag, Properties::IOFields>::initOutputModule(lowDimVtkWriter);
    lowDimProblem->addVtkOutputFields(lowDimVtkWriter);
    lowDimVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager, timeLoop, oldSol);

    // the linear solver
    using LinearSolver = BlockDiagILU0BiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // keep track of mass that left the system
    double massLeft = 0.0;

    // output file
    const auto outFileName = getParam<std::string>("Problem.OutFile");
    std::ofstream outFile(outFileName, std::ios::out);

    outFile << "[3D] hMax: " << bulkHMax << ", hMin: " << bulkHMin << std::endl;
    outFile << "[1D] hMax: " << ldHMax << ", hMin: " << ldHMin << std::endl << std::endl;

    double lowDimMass = computeGlobalMass(*lowDimProblem, sol[lowDimIdx], *lowDimGridVariables);
    double bulkMass = computeGlobalMass(*bulkProblem, sol[bulkIdx], *bulkGridVariables);
    const double initialMass = lowDimMass + bulkMass;

    std::cout << "\033[1;33m" << "The domain initially contains " << (lowDimMass + bulkMass)*1e12 << " ng tracer"
              << " (root: " << lowDimMass*1e12 << ", soil: " << bulkMass*1e12 << ")\033[0m" << '\n';

    // time loop
    timeLoop->setPeriodicCheckPoint(episodeLength);
    timeLoop->start();
    while (!timeLoop->finished())
    {
        std::cout << '\n' << "\033[1m" << "Simulation time in hours: " << timeLoop->time()/3600 << "\033[0m\n\n";

        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        oldSol = sol;
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        lowDimMass = computeGlobalMass(*lowDimProblem, sol[lowDimIdx], *lowDimGridVariables);
        bulkMass = computeGlobalMass(*bulkProblem, sol[bulkIdx], *bulkGridVariables);
        massLeft += computeGlobalBoundaryMass(*lowDimProblem, sol[lowDimIdx], *lowDimGridVariables, timeLoop->timeStepSize());

        std::cout << "\033[1;33m" << "The domain contains " << (lowDimMass + bulkMass)*1e12 << " ng tracer"
                  << " (root: " << lowDimMass*1e12 << ", soil: " << bulkMass*1e12 << ")\033[0m" << '\n';

        std::cout << "\033[1;33m" << massLeft*1e12 << " ng left domain over the root collar -> "
                  << ((lowDimMass + bulkMass) + massLeft)*1e12 << " ng balanced.\033[0m" << '\n';

        std::cout << "\033[1;33m" << "Global mass balance error: "
                  << (lowDimMass + bulkMass + massLeft - initialMass)*1e12 << " ng.\033[0m" << '\n';

        outFile << timeLoop->time() << std::scientific << std::setprecision(8)
                << " " << massLeft*1e12 << " " << lowDimMass + bulkMass + massLeft - initialMass << '\n';

        // write vtk output
        if (outputVtk && (timeLoop->isCheckPoint() || timeLoop->finished()))
        {
            bulkVtkWriter.write(timeLoop->time());
            lowDimVtkWriter.write(timeLoop->time());
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }

    outFile.close();

    timeLoop->finalize(mpiHelper.getCommunication());

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
