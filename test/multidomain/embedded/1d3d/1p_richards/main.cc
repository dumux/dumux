// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the 1d-3d embedded mixed-dimension model coupling two
 *        one-phase porous medium flow problems.
 */

#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "properties.hh"

template<class CouplingManager, class Writer>
void addElementMarkerForProjection(const CouplingManager& couplingManager, Writer& w)
{
    if constexpr(CouplingManager::couplingMode == Dumux::Embedded1d3dCouplingMode::projection)
        w.addField(couplingManager.bulkElementMarker(), "element marker");
}

// The coupling area is the area of the coupled bulk boundary faces. Each coupled face carries
// quadrature weights that sum to the face area, and a second initialization does not change anything.
template<class CouplingManager, class BulkGridGeometry>
void checkProjectionPointSources(CouplingManager& couplingManager, const BulkGridGeometry& bulkGridGeometry)
{
    if constexpr(CouplingManager::couplingMode == Dumux::Embedded1d3dCouplingMode::projection)
    {
        using namespace Dumux;
        const auto couplingArea = [&]{
            double area = 0.0;
            for (const auto& s : couplingManager.bulkPointSources())
                area += s.quadratureWeight()*s.integrationElement();
            return area;
        };

        // Per element, the weights of its point sources sum to the area of its coupled boundary faces.
        // A boundary face is coupled if a point source lies on its plane and on no other boundary face plane
        // of the element (quadrature points may also lie on an edge shared by two faces).
        const auto coupledFaceArea = [&]{
            std::vector<std::vector<const std::decay_t<decltype(couplingManager.bulkPointSources()[0])>*>> sourcesInElement(bulkGridGeometry.gridView().size(0));
            for (const auto& s : couplingManager.bulkPointSources())
                sourcesInElement[s.elementIndex()].push_back(&s);

            double totalArea = 0.0;
            for (const auto& element : elements(bulkGridGeometry.gridView()))
            {
                const auto eIdx = bulkGridGeometry.elementMapper().index(element);
                if (sourcesInElement[eIdx].empty())
                    continue;

                double weight = 0.0, area = 0.0;
                for (const auto* s : sourcesInElement[eIdx])
                    weight += s->quadratureWeight()*s->integrationElement();

                for (const auto& intersection : intersections(bulkGridGeometry.gridView(), element))
                {
                    if (!intersection.boundary())
                        continue;

                    const auto faceArea = intersection.geometry().volume();
                    const auto isOnPlane = [&](const auto& s, const auto& is) {
                        return std::abs((s.position() - is.geometry().center())*is.centerUnitOuterNormal()) < 1e-9*std::sqrt(faceArea);
                    };

                    const bool isCoupledFace = std::any_of(sourcesInElement[eIdx].begin(), sourcesInElement[eIdx].end(), [&](const auto* s) {
                        if (!isOnPlane(*s, intersection))
                            return false;
                        for (const auto& other : intersections(bulkGridGeometry.gridView(), element))
                            if (other.boundary() && other.indexInInside() != intersection.indexInInside() && isOnPlane(*s, other))
                                return false;
                        return true;
                    });

                    if (isCoupledFace)
                        area += faceArea;
                }

                if (std::abs(weight - area) > 1e-12*area)
                    DUNE_THROW(Dune::InvalidStateException, "The point source weights of element " << eIdx << " sum to "
                               << weight << ", the area of its coupled faces is " << area);
                totalArea += area;
            }
            return totalArea;
        };

        const auto numSources = couplingManager.bulkPointSources().size();
        const auto area = couplingArea();
        if (numSources == 0 || std::abs(area - coupledFaceArea()) > 1e-12*area)
            DUNE_THROW(Dune::InvalidStateException, "Coupling area " << area << " is not the area of the coupled faces " << coupledFaceArea());

        const auto localAreaFactor = couplingManager.localAreaFactor();
        couplingManager.computePointSourceData();

        if (couplingManager.bulkPointSources().size() != numSources
            || couplingManager.lowDimPointSources().size() != numSources
            || std::as_const(couplingManager).pointSourceData().size() != numSources)
            DUNE_THROW(Dune::InvalidStateException, "Number of point sources changed on re-initialization from " << numSources
                       << " to " << couplingManager.bulkPointSources().size());

        if (std::abs(couplingArea() - area) > 1e-12*area)
            DUNE_THROW(Dune::InvalidStateException, "Coupling area changed on re-initialization from " << area << " to " << couplingArea());

        const auto& newLocalAreaFactor = couplingManager.localAreaFactor();
        for (std::size_t i = 0; i < localAreaFactor.size(); ++i)
            if (std::abs(newLocalAreaFactor[i] - localAreaFactor[i]) > 1e-12)
                DUNE_THROW(Dune::InvalidStateException, "Local area factor of segment " << i << " changed on re-initialization");

        for (std::size_t id = 0; id < numSources; ++id)
            if (!(couplingManager.averageDistance(id) >= 0.0))
                DUNE_THROW(Dune::InvalidStateException, "Invalid average distance for point source " << id);

        std::cout << "Checked " << numSources << " projection point sources, coupling area " << area << std::endl;
    }
}

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
    using BulkTypeTag = Properties::TTag::SOILTYPETAG;
    using LowDimTypeTag = Properties::TTag::Root;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using BulkGridManager = Dumux::GridManager<GetPropType<BulkTypeTag, Properties::Grid>>;
    BulkGridManager bulkGridManager;
    bulkGridManager.init("Soil"); // pass parameter group

    using LowDimGridManager = Dumux::GridManager<GetPropType<LowDimTypeTag, Properties::Grid>>;
    LowDimGridManager lowDimGridManager;
    lowDimGridManager.init("Root"); // pass parameter group

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

    // the low dim spatial parameters
    using LowDimSpatialParams = GetPropType<LowDimTypeTag, Properties::SpatialParams>;
    auto lowDimSpatialParams = std::make_shared<LowDimSpatialParams>(lowDimFvGridGeometry, lowDimGridManager.getGridData());

    // the low dim problem (initial and boundary conditions)
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkFvGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(sol[bulkIdx]);
    lowDimProblem->applyInitialSolution(sol[lowDimIdx]);
    auto oldSol = sol;

    couplingManager->init(bulkProblem, lowDimProblem, sol);
    checkProjectionPointSources(*couplingManager, *bulkFvGridGeometry);
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

    // initialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(bulkVtkWriter);
    addElementMarkerForProjection(*couplingManager, bulkVtkWriter);
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

    // time loop
    timeLoop->setPeriodicCheckPoint(episodeLength);
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        oldSol = sol;
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance the time loop to the next step
        timeLoop->advanceTimeStep();

        // output the source terms
        bulkProblem->computeSourceIntegral(sol[bulkIdx], *bulkGridVariables);
        lowDimProblem->computeSourceIntegral(sol[lowDimIdx], *lowDimGridVariables);

        // write vtk output
        if (timeLoop->isCheckPoint() || timeLoop->finished())
        {
            bulkVtkWriter.write(timeLoop->time());
            lowDimVtkWriter.write(timeLoop->time());
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }

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
