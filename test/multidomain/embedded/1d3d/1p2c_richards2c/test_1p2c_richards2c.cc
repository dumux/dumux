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
 * \brief Test for the 1d-3d embedded mixed-dimension model coupling two
 *        one-phase porous medium flow problems
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
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/mixeddimension/embedded/cellcentered/bboxtreecouplingmanager.hh>
#include <dumux/mixeddimension/embedded/cellcentered/bboxtreecouplingmanagersimple.hh>
#include <dumux/mixeddimension/embedded/integrationpointsource.hh>

#include "rootproblem.hh"
#include "soilproblem.hh"

namespace Dumux {
namespace Properties {

SET_PROP(SoilTypeTag, CouplingManager)
{
    using Traits = MultiDomainTraits<TypeTag, TTAG(RootTypeTag)>;
    using type = Dumux::CCBBoxTreeEmbeddedCouplingManager<Traits>;
    // using type = Dumux::CCBBoxTreeEmbeddedCouplingManagerSimple<Traits>;
};

SET_PROP(RootTypeTag, CouplingManager)
{
    using Traits = MultiDomainTraits<TTAG(SoilTypeTag), TypeTag>;
    using type = Dumux::CCBBoxTreeEmbeddedCouplingManager<Traits>;
    // using type = Dumux::CCBBoxTreeEmbeddedCouplingManagerSimple<Traits>;
};

SET_TYPE_PROP(SoilTypeTag, PointSource, IntegrationPointSource<TypeTag>);
SET_TYPE_PROP(RootTypeTag, PointSource, IntegrationPointSource<TypeTag>);

SET_TYPE_PROP(SoilTypeTag, PointSourceHelper, IntegrationPointSourceHelper);
SET_TYPE_PROP(RootTypeTag, PointSourceHelper, IntegrationPointSourceHelper);

} // end namespace Properties

//! helper function for mass balance evaluations
template<class Problem, class SolutionVector, class GridVariables>
double computeSourceIntegral(const Problem& problem, const SolutionVector& sol, const GridVariables& gridVars)
{
    const auto& gg = problem.fvGridGeometry();
    typename SolutionVector::block_type source(0.0);
    for (const auto& element : elements(gg.gridView()))
    {
        auto fvGeometry = localView(gg);
        fvGeometry.bindElement(element);

        auto elemVolVars = localView(gridVars.curGridVolVars());
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (auto&& scv : scvs(fvGeometry))
        {
            auto pointSources = problem.scvPointSources(element, fvGeometry, elemVolVars, scv);
            // conversion to kg/s
            const auto& volVars = elemVolVars[scv];
            pointSources *= scv.volume()*volVars.extrusionFactor()
                            * volVars.density(Problem::Indices::wPhaseIdx) / volVars.molarDensity(Problem::Indices::wPhaseIdx);

            source += pointSources;
        }
    }

    std::cout << "Global integrated source (" << problem.name() << "): " << source[Problem::Indices::conti0EqIdx] << " (kg/s) / "
              <<                           source[Problem::Indices::conti0EqIdx]*3600*24*1000 << " (g/day)" << '\n';

    return source[Problem::Indices::conti0EqIdx];
}

//! helper function for mass balance evaluations
template<class Problem, class SolutionVector, class GridVariables>
double computeGlobalMass(const Problem& problem, const SolutionVector& sol, const GridVariables& gridVars)
{
    static constexpr int wPhaseIdx = Problem::Indices::wPhaseIdx;
    static constexpr int transportCompIdx = Problem::Indices::transportCompIdx;
    double mass = 0.0;

    const auto& gg = problem.fvGridGeometry();
    for (const auto& element : elements(gg.gridView()))
    {
        auto fvGeometry = localView(gg);
        fvGeometry.bindElement(element);

        auto elemVolVars = localView(gridVars.curGridVolVars());
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            mass += volVars.massFraction(wPhaseIdx, transportCompIdx)*volVars.density(wPhaseIdx)
                     *scv.volume() * volVars.porosity() * volVars.saturation(wPhaseIdx) * volVars.extrusionFactor();
        }
    }

    return mass;
}

//! helper function for mass balance evaluations
template<class Problem, class SolutionVector, class GridVariables>
double computeGlobalBoundaryMass(const Problem& problem, const SolutionVector& sol, const GridVariables& gridVars, double dt)
{
    static constexpr int transportCompIdx = Problem::Indices::transportCompIdx;
    static constexpr int transportEqIdx = Problem::Indices::transportEqIdx;
    double mass = 0.0;

    const auto& gg = problem.fvGridGeometry();
    for (const auto& element : elements(gg.gridView()))
    {
        auto fvGeometry = localView(gg);
        fvGeometry.bindElement(element);

        auto elemVolVars = localView(gridVars.curGridVolVars());
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (auto&& scvf : scvfs(fvGeometry))
            if (scvf.boundary())
                mass += problem.neumann(element, fvGeometry, elemVolVars, scvf)[transportEqIdx]
                        * scvf.area() * elemVolVars[scvf.insideScvIdx()].extrusionFactor()
                        * Problem::FluidSystem::molarMass(transportCompIdx)
                        * dt;
    }

    return mass;
}

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
    using BulkTypeTag = TTAG(SoilTypeTag);
    using LowDimTypeTag = TTAG(RootTypeTag);

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using BulkGridCreator = typename GET_PROP_TYPE(BulkTypeTag, GridCreator);
    BulkGridCreator::makeGrid("Soil"); // pass parameter group

    using LowDimGridCreator = typename GET_PROP_TYPE(LowDimTypeTag, GridCreator);
    LowDimGridCreator::makeGrid("Root"); // pass parameter group

    // we compute on the leaf grid view
    const auto& bulkGridView = BulkGridCreator::grid().leafGridView();
    const auto& lowDimGridView = LowDimGridCreator::grid().leafGridView();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // create the finite volume grid geometry
    using BulkFVGridGeometry = typename GET_PROP_TYPE(BulkTypeTag, FVGridGeometry);
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    bulkFvGridGeometry->update();
    using LowDimFVGridGeometry = typename GET_PROP_TYPE(LowDimTypeTag, FVGridGeometry);
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);
    lowDimFvGridGeometry->update();

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template DomainIdx<0>();
    constexpr auto lowDimIdx = Traits::template DomainIdx<1>();

    // the coupling manager
    using CouplingManager = typename GET_PROP_TYPE(BulkTypeTag, CouplingManager);
    auto couplingManager = std::make_shared<CouplingManager>(bulkFvGridGeometry, lowDimFvGridGeometry);

    // the problem (initial and boundary conditions)
    using BulkProblem = typename GET_PROP_TYPE(BulkTypeTag, Problem);
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, couplingManager);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimTypeTag, Problem);
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, couplingManager);

    // locally refine levels deep around the embedded grid
    int levels = getParam<int>("Soil.Grid.LocalRefinement");
    for (int i = 0; i < levels; ++i)
    {
        auto& soilGrid = BulkGridCreator::grid();
        using BulkGridView = typename GET_PROP_TYPE(BulkTypeTag, GridView);
        using LowDimGridView = typename GET_PROP_TYPE(LowDimTypeTag, GridView);

        CCMixedDimensionGlue<BulkGridView, LowDimGridView>
            glue(bulkFvGridGeometry->boundingBoxTree(), lowDimFvGridGeometry->boundingBoxTree());

        // refine all 3D cells intersected
        for (const auto& is : intersections(glue))
        {
            for (unsigned int outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                const auto cutElement = is.outside(outsideIdx);

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

        // update the bounding box tree
        bulkFvGridGeometry->update();
    }

    // update geometry after refinement
    bulkFvGridGeometry->update();
    lowDimFvGridGeometry->update();

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

    std::cout << "[3D] hMax: " << bulkHMax << ", hMin: " << bulkHMin << std::endl;
    std::cout << "[1D] hMax: " << ldHMax << ", hMin: " << ldHMin << std::endl;

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
    using BulkGridVariables = typename GET_PROP_TYPE(BulkTypeTag, GridVariables);
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    bulkGridVariables->init(sol[bulkIdx], oldSol[bulkIdx]);
    using LowDimGridVariables = typename GET_PROP_TYPE(LowDimTypeTag, GridVariables);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    lowDimGridVariables->init(sol[lowDimIdx], oldSol[lowDimIdx]);

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto episodeLength = getParam<Scalar>("TimeLoop.EpisodeLength");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    const bool outputVtk = getParam<bool>("Problem.EnableVtkOutput", true);

    // intialize the vtk output module
    VtkOutputModule<BulkTypeTag> bulkVtkWriter(*bulkProblem, *bulkFvGridGeometry, *bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GET_PROP_TYPE(BulkTypeTag, VtkOutputFields)::init(bulkVtkWriter);
    if (outputVtk) bulkVtkWriter.write(0.0);

    VtkOutputModule<LowDimTypeTag> lowDimVtkWriter(*lowDimProblem, *lowDimFvGridGeometry, *lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
    GET_PROP_TYPE(LowDimTypeTag, VtkOutputFields)::init(lowDimVtkWriter);
    lowDimProblem->addVtkOutputFields(lowDimVtkWriter);
    if (outputVtk) lowDimVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager, timeLoop);

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

    const double initialMass = computeGlobalMass(*lowDimProblem, sol[lowDimIdx], *lowDimGridVariables) +
                               computeGlobalMass(*bulkProblem, sol[bulkIdx], *bulkGridVariables);

    // time loop
    timeLoop->setPeriodicCheckPoint(episodeLength);
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(oldSol);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        oldSol = sol;
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // output the source terms
        // computeSourceIntegral(*bulkProblem, sol[bulkIdx], *bulkGridVariables);
        // computeSourceIntegral(*lowDimProblem, sol[lowDimIdx], *lowDimGridVariables);

        double lowDimMass = computeGlobalMass(*lowDimProblem, sol[lowDimIdx], *lowDimGridVariables);
        double bulkMass = computeGlobalMass(*bulkProblem, sol[bulkIdx], *bulkGridVariables);
        massLeft += computeGlobalBoundaryMass(*lowDimProblem, sol[lowDimIdx], *lowDimGridVariables, timeLoop->timeStepSize());

        std::cout << "\033[1;33m" << "The domain contains " << (lowDimMass + bulkMass)*1e12 << " ng tracer"
                  << " (root: " << lowDimMass*1e12 << ", soil: " << bulkMass*1e12 << ")\033[0m" << '\n';

        std::cout << "\033[1;33m" << massLeft*1e12 << " ng left domain over the root collar -> "
                  << ((lowDimMass + bulkMass) + massLeft)*1e12 << " ng balanced.\033[0m" << '\n';

        std::cout << "\033[1;33m" << " Global mass balance error: "
                  << (lowDimMass + bulkMass + massLeft - initialMass)*1e12 << " ng.\033[0m" << '\n';

        outFile << timeLoop->time() << std::scientific << std::setprecision(8)
                << " " << massLeft*1e12 << " " << lowDimMass + bulkMass + massLeft - initialMass << '\n';

        // write vtk output
        if (outputVtk)
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

    timeLoop->finalize(mpiHelper.getCollectiveCommunication());

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
