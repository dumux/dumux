// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The block-diagonal AMG solver of the embedded coupling solves the coupled system on distributed grids
 *
 * The problem is linear, with Dirichlet conditions on the tissue boundary. Every rank solves it on the whole grids
 * with a direct solver as the reference. On the distributed grids, the residual norm of the solver has to equal the
 * norm of the sequential residual for the problem without Dirichlet conditions, whose residual is a sum of the
 * partial residuals of the processes in every row; the state after one Newton step with the iterative solver, started
 * from an initial guess that differs between the processes holding a degree of freedom, has to have a residual
 * below the tolerance relative to the initial one and has to equal the reference; the iterative solver may need at
 * most 1.5 times plus 3 the iterations it needs on the whole grids; and Newton's method with the iterative solver has
 * to converge to the reference.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/blockdiagamgsolver.hh>

#include "properties.hh"
#include "setup.hh"

namespace Dumux {

using BulkTypeTag = Properties::TTag::BULKTYPETAG;
using LowDimTypeTag = Properties::TTag::LOWDIMTYPETAG;
using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();
using TestCouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
using Assembler = MultiDomainFVAssembler<Traits, TestCouplingManager, DiffMethod::numeric>;
using SolutionVector = typename Traits::SolutionVector;
using Communication = Dune::Communication<Dune::MPIHelper::MPICommunicator>;

//! The coupled problem on two grid views, optionally with Dirichlet conditions on the tissue boundary
struct CoupledProblem
{
    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;

    template<class BulkGridView, class LowDimGridView>
    CoupledProblem(const BulkGridView& bulkGridView, const LowDimGridView& lowDimGridView, bool dirichlet = true)
    {
        bulkGridGeometry = std::make_shared<BulkGridGeometry>(bulkGridView);
        lowDimGridGeometry = std::make_shared<LowDimGridGeometry>(lowDimGridView);
        couplingManager = std::make_shared<TestCouplingManager>(bulkGridGeometry, lowDimGridGeometry);
        couplingManager->setLowDimElementData(1, [gg = lowDimGridGeometry](std::size_t eIdx, std::span<double> values)
        { values[0] = 1.0 + gg->element(eIdx).geometry().center()[1]; });
        couplingManager->setPointSourceValues(1, [](const auto& source, std::span<double> values)
        { values[0] = 1.0 + 0.5*source.quadratureWeight(); });

        bulkProblem = std::make_shared<BulkProblem>(bulkGridGeometry, couplingManager);
        bulkProblem->setDirichletBoundary(dirichlet);
        lowDimProblem = std::make_shared<LowDimProblem>(lowDimGridGeometry, couplingManager);

        sol[bulkIdx].resize(bulkGridGeometry->numDofs());
        sol[lowDimIdx].resize(lowDimGridGeometry->numDofs());
        const auto set = [](const auto& gridGeometry, auto& v, auto&& f)
        {
            auto fvGeometry = localView(gridGeometry);
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                    v[scv.dofIndex()] = f(scv.dofPosition());
            }
        };
        set(*bulkGridGeometry, sol[bulkIdx], [](const auto& x){ return bulkPressure(x, 0); });
        set(*lowDimGridGeometry, sol[lowDimIdx], [](const auto& x){ return lowDimPressure(x, 0); });

        couplingManager->init(bulkProblem, lowDimProblem, sol);
        bulkProblem->computePointSourceMap();
        lowDimProblem->computePointSourceMap();

        bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGridGeometry);
        lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGridGeometry);
        bulkGridVariables->init(sol[bulkIdx]);
        lowDimGridVariables->init(sol[lowDimIdx]);

        assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                std::make_tuple(bulkGridGeometry, lowDimGridGeometry),
                                                std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                couplingManager);
    }

    auto gridGeometries() const
    { return std::make_tuple(bulkGridGeometry, lowDimGridGeometry); }

    std::shared_ptr<BulkGridGeometry> bulkGridGeometry;
    std::shared_ptr<LowDimGridGeometry> lowDimGridGeometry;
    std::shared_ptr<TestCouplingManager> couplingManager;
    std::shared_ptr<BulkProblem> bulkProblem;
    std::shared_ptr<LowDimProblem> lowDimProblem;
    std::shared_ptr<BulkGridVariables> bulkGridVariables;
    std::shared_ptr<LowDimGridVariables> lowDimGridVariables;
    std::shared_ptr<Assembler> assembler;
    SolutionVector sol;
};

//! The dofs whose rows this process assembles and their positions
template<class GridGeometry>
std::vector<std::pair<std::size_t, DofKey>> assembledDofs(const GridGeometry& gridGeometry)
{
    std::vector<std::pair<std::size_t, DofKey>> dofs;
    const auto& gridView = gridGeometry.gridView();
    if constexpr (GridGeometry::discMethod == DiscretizationMethods::box)
    {
        for (const auto& vertex : vertices(gridView))
            if (vertex.partitionType() == Dune::InteriorEntity || vertex.partitionType() == Dune::BorderEntity)
                dofs.emplace_back(gridGeometry.vertexMapper().index(vertex), dofKey(vertex.geometry().center()));
    }
    else
    {
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
            dofs.emplace_back(gridGeometry.elementMapper().index(element), dofKey(element.geometry().center()));
    }
    return dofs;
}

//! The values of a solution by the positions of the dofs whose rows this process assembles
std::array<std::map<DofKey, double>, 2> valuesByPosition(const CoupledProblem& problem, const SolutionVector& x)
{
    std::array<std::map<DofKey, double>, 2> values;
    for (const auto& [dof, key] : assembledDofs(*problem.bulkGridGeometry))
        values[0][key] = x[bulkIdx][dof][0];
    for (const auto& [dof, key] : assembledDofs(*problem.lowDimGridGeometry))
        values[1][key] = x[lowDimIdx][dof][0];
    return values;
}

//! The largest difference to the reference over all processes, relative to the largest reference value
double relativeDifference(const std::array<std::map<DofKey, double>, 2>& values,
                          const std::array<std::map<DofKey, double>, 2>& reference,
                          const Communication& comm)
{
    double difference = 0.0, magnitude = 0.0;
    for (int domain = 0; domain < 2; ++domain)
    {
        for (const auto& [key, value] : reference[domain])
            magnitude = std::max(magnitude, std::abs(value));
        for (const auto& [key, value] : values[domain])
        {
            const auto it = reference[domain].find(key);
            if (it == reference[domain].end())
                DUNE_THROW(Dune::InvalidStateException, "A dof is not in the reference");
            if (!std::isfinite(value))
                difference = std::numeric_limits<double>::infinity();
            difference = std::max(difference, std::abs(value - it->second));
        }
    }
    return comm.max(difference)/magnitude;
}

//! The direct solver of the whole grids, which every process of a parallel run uses with its own copy of the grids
class ReferenceSolver : public UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>
{
    using ParentType = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
public:
    using ParentType::ParentType;

    template<class Vector>
    auto norm(const Vector& x) const
    { return x.two_norm(); }
};

template<class LinearSolver>
class StepCountingNewtonSolver : public MultiDomainNewtonSolver<Assembler, LinearSolver, TestCouplingManager>
{
    using ParentType = MultiDomainNewtonSolver<Assembler, LinearSolver, TestCouplingManager>;
public:
    using ParentType::ParentType;

    int numSteps() const
    { return this->numSteps_; }
};

template<class LinearSolver>
int runTest(const Communication& comm, const Communication& self)
{
    using BulkGrid = GetPropType<BulkTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimTypeTag, Properties::Grid>;
    const auto& mpiHelper = Dune::MPIHelper::instance();

    int failures = 0;
    const auto fail = [&](const auto& message){ std::cerr << "rank " << comm.rank() << ": " << message << std::endl; ++failures; };

    // the whole grids with a sequential communicator: the iterations of the iterative solver and the reference
    const auto wholeBulkGrid = makeBulkGrid<BulkGrid>(mpiHelper.getLocalCommunicator());
    const auto wholeLowDimGrid = makeNetwork<LowDimGrid>();
    CoupledProblem whole(wholeBulkGrid->leafGridView(), wholeLowDimGrid->leafGridView());
    whole.assembler->assembleJacobianAndResidual(whole.sol);
    LinearSolver wholeSolver(whole.gridGeometries(), whole.couplingManager);
    auto wholeUpdate = whole.sol;
    wholeUpdate = 0.0;
    const auto wholeResult = wholeSolver.solve(whole.assembler->jacobian(), wholeUpdate, whole.assembler->residual());

    auto referenceSolution = whole.sol;
    StepCountingNewtonSolver<ReferenceSolver> wholeNewton(whole.assembler, std::make_shared<ReferenceSolver>(), whole.couplingManager, self);
    wholeNewton.solve(referenceSolution);
    const auto reference = valuesByPosition(whole, referenceSolution);

    // one Newton step on the distributed grids
    const auto bulkGrid = makeBulkGrid<BulkGrid>(mpiHelper.getCommunicator());
    auto lowDimGrid = makeNetwork<LowDimGrid>();
    distributeNetwork(*lowDimGrid, bulkGrid->leafGridView());
    CoupledProblem distributed(bulkGrid->leafGridView(), lowDimGrid->leafGridView());
    distributed.assembler->assembleJacobianAndResidual(distributed.sol);
    auto linearSolver = std::make_shared<LinearSolver>(distributed.gridGeometries(), distributed.couplingManager);
    const auto initialNorm = linearSolver->norm(distributed.assembler->residual());

    CoupledProblem wholeNeumann(wholeBulkGrid->leafGridView(), wholeLowDimGrid->leafGridView(), false);
    CoupledProblem distributedNeumann(bulkGrid->leafGridView(), lowDimGrid->leafGridView(), false);
    wholeNeumann.assembler->assembleJacobianAndResidual(wholeNeumann.sol);
    distributedNeumann.assembler->assembleJacobianAndResidual(distributedNeumann.sol);
    const auto wholeNorm = wholeNeumann.assembler->residual().two_norm();
    const auto distributedNorm = linearSolver->norm(distributedNeumann.assembler->residual());
    auto update = distributed.sol;
    update = 1.0 + comm.rank();
    const auto result = linearSolver->solve(distributed.assembler->jacobian(), update, distributed.assembler->residual());

    auto state = distributed.sol;
    state -= update;
    distributed.couplingManager->updateSolution(state);
    auto residual = distributed.assembler->residual();
    distributed.assembler->assembleResidual(residual, state);
    const auto residualReduction = linearSolver->norm(residual)/initialNorm;
    const auto stepDifference = relativeDifference(valuesByPosition(distributed, state), reference, comm);

    // Newton's method with the iterative solver
    CoupledProblem distributedNewton(bulkGrid->leafGridView(), lowDimGrid->leafGridView());
    auto newtonSolution = distributedNewton.sol;
    StepCountingNewtonSolver<LinearSolver> newton(distributedNewton.assembler,
        std::make_shared<LinearSolver>(distributedNewton.gridGeometries(), distributedNewton.couplingManager),
        distributedNewton.couplingManager, comm);
    int newtonSteps = -1;
    try
    {
        newton.solve(newtonSolution);
        newtonSteps = newton.numSteps();
    }
    catch (const NumericalProblem&) {}
    const auto newtonDifference = relativeDifference(valuesByPosition(distributedNewton, newtonSolution), reference, comm);

    if (comm.rank() == 0)
        std::cout << "residual norm without Dirichlet conditions: distributed " << distributedNorm << ", sequential " << wholeNorm << "\n"
                  << "iterations: " << result.iterations << " distributed, " << wholeResult.iterations << " on the whole grids\n"
                  << "one step: |residual| / |initial residual| = " << residualReduction
                  << ", max |state - reference| / max |reference| = " << stepDifference << "\n"
                  << "Newton: " << newtonSteps << " steps (-1: not converged), max |solution - reference| / max |reference| = "
                  << newtonDifference << std::endl;

    if (std::abs(distributedNorm - wholeNorm) > 1e-12*wholeNorm)
        fail("the residual norm differs from the norm of the sequential residual");
    if (!result.converged || !wholeResult.converged)
        fail("the iterative solver did not converge");
    if (result.iterations > 1.5*wholeResult.iterations + 3)
        fail("the iterative solver needs more iterations on the distributed grids than on the whole grids");
    if (residualReduction > 1e-7)
        fail("the residual after one step is not reduced to the tolerance");
    if (stepDifference > 1e-6)
        fail("the state after one step differs from the reference");
    if (newtonSteps < 0)
        fail("Newton's method with the iterative solver did not converge");
    if (newtonDifference > 1e-9)
        fail("the solution of Newton's method differs from the reference");

    return failures;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    const Communication comm(mpiHelper.getCommunicator());
    const Communication self(mpiHelper.getLocalCommunicator());

    const auto krylov = getParam<std::string>("Test.Krylov", "BiCGSTAB");
    int failures = 0;
    if (krylov == "BiCGSTAB")
        failures = runTest<EmbeddedBlockDiagAMGBiCGSTABSolver<TestCouplingManager>>(comm, self);
    else if (krylov == "GMRes")
        failures = runTest<EmbeddedBlockDiagAMGRestartedGMResSolver<TestCouplingManager>>(comm, self);
    else
        DUNE_THROW(Dune::InvalidStateException, "Unknown Krylov solver " << krylov);

    failures = comm.sum(failures);
    if (comm.rank() == 0)
        std::cout << (failures == 0 ? "PASSED" : "FAILED") << " on " << comm.size() << " processes (" << krylov << ")" << std::endl;

    return failures == 0 ? 0 : 1;
}
