// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Test for the 1d-3d embedded mixed-dimension model coupling two one-phase
 *        problems stated with the variables of the local dofs, i.e. also for hybrid
 *        and finite element schemes.
 */

#include <config.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/geometry/diameter.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/parallel/multithreading.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/embedded/elementdofs.hh>

#include <test/testing/convergencetest.hh>
#include <test/multidomain/embedded/jacobiancheck.hh>
#include <test/multidomain/embedded/coloringcheck.hh>

#include "properties.hh"

//! The relative L2 error of a solution with the basis functions of the grid discretization
template<class GridDiscretization, class SolutionVector, class ExactSolution>
double relativeL2Error(const GridDiscretization& gridDiscretization, const SolutionVector& sol,
                       ExactSolution&& exactSolution, int order)
{
    using namespace Dumux;
    double error = 0.0, volume = 0.0;
    auto elemDisc = localView(gridDiscretization);
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        elemDisc.bindElement(element);
        const auto dofIndices = EmbeddedCoupling::elementDofIndices(elemDisc);
        const auto geometry = element.geometry();
        const auto& quad = Dune::QuadratureRules<double, GridDiscretization::GridView::dimension>::rule(geometry.type(), order);
        for (const auto& qp : quad)
        {
            const auto globalPos = geometry.global(qp.position());
            std::vector<Dune::FieldVector<double, 1>> shapeValues;
            EmbeddedCoupling::basisFunctionValues(gridDiscretization, geometry, globalPos, shapeValues);
            double value = 0.0;
            for (std::size_t i = 0; i < dofIndices.size(); ++i)
                value += sol[dofIndices[i]][0]*shapeValues[i][0];
            const auto weight = qp.weight()*geometry.integrationElement(qp.position());
            error += (value - exactSolution(globalPos))*(value - exactSolution(globalPos))*weight;
            volume += weight;
        }
    }
    return std::sqrt(error)/volume;
}

//! The sum of the values of the point sources, each point counted once
template<class CouplingManager, std::size_t i, class Problem, class SolutionVector>
double totalPointSourceValue(const CouplingManager& couplingManager, Dune::index_constant<i> domainIdx,
                             const Problem& problem, const SolutionVector& sol)
{
    using namespace Dumux;
    const auto& gridDiscretization = Dumux::gridDiscretization(problem);
    auto elemDisc = localView(gridDiscretization);
    const auto& pointSources = problem.pointSources();
    double total = 0.0;
    for (const auto& context : couplingManager.pointSourceContexts(domainIdx))
    {
        elemDisc.bindElement(gridDiscretization.element(context.elementIndex()));
        total += pointSources.eval(elemDisc, 0, context)[0];
    }
    return total;
}

template<class GridDiscretization>
constexpr bool hasSubControlVolumes = requires (const typename GridDiscretization::LocalView& elemDisc) { scvs(elemDisc); };

//! The sum of the values of the point sources that act on the control volumes
template<class Problem>
double totalControlVolumeSource(const Problem& problem)
{
    using namespace Dumux;
    const auto& gridDiscretization = Dumux::gridDiscretization(problem);
    auto elemDisc = localView(gridDiscretization);
    const auto& pointSources = problem.pointSources();
    double total = 0.0;
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        elemDisc.bindElement(element);
        if constexpr (hasSubControlVolumes<std::decay_t<decltype(gridDiscretization)>>)
            for (const auto& scv : scvs(elemDisc))
                for (const auto& context : pointSources.contexts(elemDisc, scv))
                    total += pointSources.eval(elemDisc, 0, context)[0];
    }
    return total;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    using BulkTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimTypeTag = Properties::TTag::LOWDIMTYPETAG;

    using BulkGridManager = Dumux::GridManager<GetPropType<BulkTypeTag, Properties::Grid>>;
    BulkGridManager bulkGridManager;
    bulkGridManager.init("Tissue");

    using LowDimGridManager = Dumux::GridManager<GetPropType<LowDimTypeTag, Properties::Grid>>;
    LowDimGridManager lowDimGridManager;
    lowDimGridManager.init("Vessel");

    const auto& bulkGridView = bulkGridManager.grid().leafGridView();
    const auto& lowDimGridView = lowDimGridManager.grid().leafGridView();

    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    auto bulkGridGeometry = std::make_shared<BulkGridGeometry>(bulkGridView);
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    auto lowDimGridGeometry = std::make_shared<LowDimGridGeometry>(lowDimGridView);

    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();

    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(bulkGridGeometry, lowDimGridGeometry);

    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkGridGeometry, couplingManager);
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimGridGeometry, couplingManager);

    Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimGridGeometry->numDofs());
    sol = 0.0;

    couplingManager->init(bulkProblem, lowDimProblem, sol);
    bulkProblem->setPointSources();
    lowDimProblem->setPointSources();

    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGridGeometry);
    bulkGridVariables->init(sol[bulkIdx]);
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGridGeometry);
    lowDimGridVariables->init(sol[lowDimIdx]);

    using Assembler = Experimental::MultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    const auto makeAssembler = [&]()
    {
        return std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                           std::make_tuple(bulkGridGeometry, lowDimGridGeometry),
                                           std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                           couplingManager);
    };
    auto assembler = makeAssembler();

    if (getParam<bool>("Test.CheckJacobian", false))
    {
        const auto result = Test::compareJacobianWithFiniteDifferences(
            *assembler, *couplingManager, sol, {1.0, 1e-2, 1e-4, 1e-6, 1e-8}
        );
        double best = std::numeric_limits<double>::max();
        for (const auto& [epsilon, relativeError] : result)
        {
            std::cout << "[jacobian] eps = " << epsilon << ", |J v - (R(u + eps v) - R(u))/eps|/|J v| = " << relativeError << std::endl;
            best = std::min(best, relativeError);
        }
        if (!(best < getParam<double>("Test.JacobianTolerance", 1e-6)))
            DUNE_THROW(Dune::InvalidStateException, "Jacobian does not match finite differences of the residual, best relative error " << best);
    }

    if (getParam<bool>("Test.CheckColoring", false))
    {
        couplingManager->computeColorsForAssembly();
        const auto numConflicts = Test::countColoringConflicts(*couplingManager);
        std::cout << "[coloring] conflicts between elements of the same color: " << numConflicts << std::endl;
        if (numConflicts > 0)
            DUNE_THROW(Dune::InvalidStateException, "Elements of the same color conflict");
    }

    if (getParam<bool>("Test.CheckMultithreading", false))
    {
        const auto* numThreads = std::getenv("DUMUX_NUM_THREADS");
        if (Multithreading::isSerial() || !numThreads || std::atoi(numThreads) < 2)
            DUNE_THROW(Dune::InvalidStateException, "The multithreading check needs a threading backend and DUMUX_NUM_THREADS > 1");

        // the assemblers read the parameter when they are constructed
        Parameters::init([](Dune::ParameterTree& params) { params["Assembly.Multithreading"] = "false"; });
        auto serialAssembler = makeAssembler();
        Parameters::init([](Dune::ParameterTree& params) { params["Assembly.Multithreading"] = "true"; });
        auto threadedAssembler = makeAssembler();
        couplingManager->updateSolution(sol);
        serialAssembler->assembleJacobianAndResidual(sol);
        threadedAssembler->assembleJacobianAndResidual(sol);
        const auto difference = Test::compareMatrices(serialAssembler->jacobian(), threadedAssembler->jacobian());
        auto residualDifference = serialAssembler->residual();
        residualDifference -= threadedAssembler->residual();
        std::cout << "[multithreading] max |J_serial - J_threaded| = " << difference.maxAbsDifference
                  << ", max |J| = " << difference.maxAbsEntry
                  << ", identical bit by bit: " << std::boolalpha << difference.identical
                  << ", |R_serial - R_threaded| = " << Test::norm(residualDifference) << std::endl;
        const auto tolerance = getParam<double>("Test.MultithreadingTolerance", 1e-12);
        if (difference.maxAbsDifference > tolerance*difference.maxAbsEntry
            || Test::norm(residualDifference) > tolerance*Test::norm(serialAssembler->residual()))
            DUNE_THROW(Dune::InvalidStateException, "Multithreaded and serial assembly differ");
    }

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
    nonLinearSolver.solve(sol);

    // the exchange between the domains is conserved
    const auto bulkSource = totalPointSourceValue(*couplingManager, bulkIdx, *bulkProblem, sol[bulkIdx]);
    const auto lowDimSource = totalPointSourceValue(*couplingManager, lowDimIdx, *lowDimProblem, sol[lowDimIdx]);
    const auto exchangeError = std::abs(bulkSource + lowDimSource)/std::abs(lowDimSource);
    std::cout << "[conservation] total bulk source = " << bulkSource << ", total low-dim source = " << lowDimSource
              << ", rel. error = " << exchangeError << std::endl;
    if (exchangeError > 1e-10)
        DUNE_THROW(Dune::InvalidStateException, "Exchange between the domains is not conservative, rel. error = " << exchangeError);

    // the control volumes receive the total source
    if constexpr (hasSubControlVolumes<BulkGridGeometry>)
    {
        const auto bulkCVSource = totalControlVolumeSource(*bulkProblem);
        const auto cvError = std::abs(bulkCVSource - bulkSource)/std::abs(bulkSource);
        std::cout << "[conservation] bulk source of the control volumes = " << bulkCVSource << ", rel. error = " << cvError << std::endl;
        if (cvError > 1e-10)
            DUNE_THROW(Dune::InvalidStateException, "The control volumes do not receive the bulk source, rel. error = " << cvError);
    }
    if constexpr (hasSubControlVolumes<LowDimGridGeometry>)
    {
        const auto lowDimCVSource = totalControlVolumeSource(*lowDimProblem);
        const auto cvError = std::abs(lowDimCVSource - lowDimSource)/std::abs(lowDimSource);
        std::cout << "[conservation] low-dim source of the control volumes = " << lowDimCVSource << ", rel. error = " << cvError << std::endl;
        if (cvError > 1e-10)
            DUNE_THROW(Dune::InvalidStateException, "The control volumes do not receive the low-dim source, rel. error = " << cvError);
    }

    static const int order = getParam<int>("L2Error.IntegrationOrder");
    const auto norm3D = relativeL2Error(*bulkGridGeometry, sol[bulkIdx], [&](const auto& pos) { return bulkProblem->exactSolution(pos); }, order);
    const auto norm1D = relativeL2Error(*lowDimGridGeometry, sol[lowDimIdx], [&](const auto& pos) { return lowDimProblem->exactSolution(pos); }, order);

    double hMaxBulk = 0.0;
    for (const auto& element : elements(bulkGridView))
        hMaxBulk = std::max(hMaxBulk, diameter(element.geometry()));

    std::cout << "-----------------------------------------------------\n"
              << " L2_p1d: " << norm1D << '\n'
              << " L2_p3d: " << norm3D << " hmax_3d: " << hMaxBulk << '\n'
              << "-----------------------------------------------------" << std::endl;

    Test::writeConvergenceTestErrors(getParam<std::string>("Problem.Name"), hMaxBulk, {
        {"L2Rel(p_3d)", norm3D}, {"L2Rel(p_1d)", norm1D}
    });

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
