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
#include <random>
#include <type_traits>
#include <vector>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/referenceelements.hh>
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
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <test/testing/convergencetest.hh>

#include "properties.hh"
#include "l2norm.hh"

/*!
 * \brief Integrate the coupling kernel against unity.
 *
 *        Summing every source weight over all bulk elements (and, for box, all sub-control
 *        volumes) must reproduce the total vessel length, because the kernel is normalised to
 *        integrate to one per unit length. The check is independent of the physics and of the
 *        discretisation, so it isolates the distribution itself: if the weights do not sum
 *        correctly the coupling is losing or inventing exchange, however plausible the
 *        resulting pressure field looks.
 *
 *        Templated so that the branch is discarded for coupling modes that carry no kernel
 *        weights at all -- `if constexpr` in a non-template context would still require the
 *        member to exist.
 */
template<class CouplingManager, class BulkGG, class LowDimGG>
void checkKernelConservation(const CouplingManager& couplingManager,
                             const BulkGG& bulkGG, const LowDimGG& lowDimGG)
{
    using namespace Dumux;
    if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
    {
        constexpr int nScv = CouplingManager::maxNumBulkElementDofs();
        double total = 0.0;
        for (std::size_t eIdx = 0; eIdx < bulkGG.gridView().size(0); ++eIdx)
            for (int s = 0; s < nScv; ++s)
                for (const auto w : couplingManager.bulkSourceWeights(eIdx, s))
                    total += w;

        double length = 0.0;
        for (const auto& element : elements(lowDimGG.gridView()))
            length += element.geometry().volume();

        const auto relErr = std::abs(total - length)/length;
        std::cout << "[conservation] integral of the kernel against unity = " << total
                  << ", vessel length = " << length
                  << ", rel. error = " << relErr << std::endl;
        if (relErr > 1e-10)
            DUNE_THROW(Dune::InvalidStateException,
                       "Kernel distribution is not conservative, rel. error = " << relErr);
    }
}

/*!
 * \brief Check that the kernel mass of each bulk element is the mass of the kernel inside it.
 *
 *        The kernel is constant, 1/(pi rho^2), inside the cylinder of radius rho around each
 *        segment, so the kernel mass of a bulk element is the volume of its intersection with
 *        the cylinders divided by pi rho^2. This volume is estimated by uniform random sampling
 *        of the element, independently of the coupling manager. The check fails if the
 *        manager assigns kernel mass to elements that do not contain it.
 */
template<class CouplingManager, class BulkGG, class LowDimGG, class LowDimSpatialParams>
void checkKernelPlacement(const CouplingManager& couplingManager,
                          const BulkGG& bulkGG, const LowDimGG& lowDimGG,
                          const LowDimSpatialParams& lowDimSpatialParams)
{
    using namespace Dumux;
    if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
    {
        constexpr int dim = BulkGG::GridView::dimension;
        constexpr int nScv = CouplingManager::maxNumBulkElementDofs();
        const auto kernelWidthFactor = getParam<double>("MixedDimension.KernelWidthFactor");

        using GlobalPosition = typename BulkGG::GridView::template Codim<0>::Geometry::GlobalCoordinate;
        struct Segment { GlobalPosition a, b; double rho; };
        std::vector<Segment> segments;
        for (const auto& element : elements(lowDimGG.gridView()))
        {
            const auto geometry = element.geometry();
            const auto rho = kernelWidthFactor*lowDimSpatialParams.radius(lowDimGG.elementMapper().index(element));
            segments.push_back({geometry.corner(0), geometry.corner(1), rho});
        }

        const auto distanceToAxis = [](const GlobalPosition& p, const Segment& s, double& t)
        {
            const auto ab = s.b - s.a;
            t = (p - s.a)*ab/ab.two_norm2();
            auto proj = s.a; proj.axpy(t, ab);
            return (p - proj).two_norm();
        };

        std::mt19937 generator(42);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        constexpr int numSamples = 4000;

        double misplaced = 0.0, total = 0.0;
        for (const auto& element : elements(bulkGG.gridView()))
        {
            const auto eIdx = bulkGG.elementMapper().index(element);
            double weight = 0.0;
            for (int s = 0; s < nScv; ++s)
                for (const auto w : couplingManager.bulkSourceWeights(eIdx, s))
                    weight += w;

            const auto geometry = element.geometry();
            const auto center = geometry.center();
            double diameter = 0.0;
            for (int i = 0; i < geometry.corners(); ++i)
                diameter = std::max(diameter, 2.0*(geometry.corner(i) - center).two_norm());

            std::vector<const Segment*> nearSegments;
            for (const auto& s : segments)
            {
                double t;
                const auto length = (s.b - s.a).two_norm();
                if (distanceToAxis(center, s, t) < s.rho + diameter && t > -diameter/length && t < 1.0 + diameter/length)
                    nearSegments.push_back(&s);
            }

            double reference = 0.0;
            if (!nearSegments.empty())
            {
                const bool isSimplex = geometry.type().isSimplex();
                const auto refVolume = Dune::referenceElement(geometry).volume();
                for (int i = 0; i < numSamples; ++i)
                {
                    typename std::decay_t<decltype(geometry)>::LocalCoordinate local;
                    for (int d = 0; d < dim; ++d)
                        local[d] = uniform(generator);
                    if (isSimplex)
                    {
                        std::sort(local.begin(), local.end());
                        for (int d = dim-1; d > 0; --d)
                            local[d] -= local[d-1];
                    }

                    const auto global = geometry.global(local);
                    for (const auto* s : nearSegments)
                    {
                        double t;
                        if (distanceToAxis(global, *s, t) <= s->rho && t >= 0.0 && t <= 1.0)
                            reference += geometry.integrationElement(local)*refVolume/(M_PI*s->rho*s->rho);
                    }
                }
                reference /= numSamples;
            }

            misplaced += std::abs(weight - reference);
            total += reference;
        }

        const auto misplacedFraction = 0.5*misplaced/total;
        std::cout << "[placement] fraction of the kernel mass in wrong bulk elements = " << misplacedFraction << std::endl;
        if (misplacedFraction > 0.05)
            DUNE_THROW(Dune::InvalidStateException,
                       "Kernel mass assigned to wrong bulk elements, misplaced fraction = " << misplacedFraction);
    }
}

//! The source (point sources and distributed sources) integrated over a domain
template<class Problem, class SolutionVector, class GridVariables>
double integrateCouplingSource(const Problem& problem, const SolutionVector& sol, const GridVariables& gridVars)
{
    using namespace Dumux;
    double total = 0.0;
    auto fvGeometry = localView(problem.gridGeometry());
    auto elemVolVars = localView(gridVars.curGridVolVars());
    for (const auto& element : elements(problem.gridGeometry().gridView()))
    {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);
        for (const auto& scv : scvs(fvGeometry))
            total += (problem.scvPointSources(element, fvGeometry, elemVolVars, scv)[0]
                      + problem.source(element, fvGeometry, elemVolVars, scv)[0])
                     *scv.volume()*elemVolVars[scv].extrusionFactor();
    }
    return total;
}

//! The mass that leaves the one-dimensional domain enters the bulk domain
template<class CouplingManager, class BulkProblem, class LowDimProblem, class SolutionVector, class BulkGridVariables, class LowDimGridVariables>
void checkExchangeConservation(const BulkProblem& bulkProblem, const LowDimProblem& lowDimProblem,
                               const SolutionVector& sol,
                               const BulkGridVariables& bulkGridVariables, const LowDimGridVariables& lowDimGridVariables)
{
    using namespace Dumux;
    if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
    {
        const auto bulkSource = integrateCouplingSource(bulkProblem, sol[Dune::index_constant<0>{}], bulkGridVariables);
        const auto lowDimSource = integrateCouplingSource(lowDimProblem, sol[Dune::index_constant<1>{}], lowDimGridVariables);
        const auto relErr = std::abs(bulkSource + lowDimSource)/std::abs(lowDimSource);
        std::cout << "[conservation] total bulk source = " << bulkSource << ", total low-dim source = " << lowDimSource
                  << ", rel. error = " << relErr << std::endl;
        if (relErr > 1e-10)
            DUNE_THROW(Dune::InvalidStateException,
                       "Exchange between the domains is not conservative, rel. error = " << relErr);
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

    checkKernelConservation(*couplingManager, *bulkFvGridGeometry, *lowDimFvGridGeometry);
    checkKernelPlacement(*couplingManager, *bulkFvGridGeometry, *lowDimFvGridGeometry, lowDimProblem->spatialParams());

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

    checkExchangeConservation<CouplingManager>(*bulkProblem, *lowDimProblem, sol, *bulkGridVariables, *lowDimGridVariables);

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
    static const int order = getParam<int>("L2Error.IntegrationOrder");
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

    Test::writeConvergenceTestErrors(getParam<std::string>("Problem.Name"), hMaxBulk, {
        {"L2Rel(p_3d)", norm3D}, {"L2Rel(p_1d)", norm1D}, {"L2Rel(q)", sourceNorm}
    });

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
