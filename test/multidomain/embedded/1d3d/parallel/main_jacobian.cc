// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The Jacobian of the embedded coupling on distributed grids equals the sequential one, including the
 *        entries whose column is a degree of freedom of another process
 *
 * Every process assembles the Jacobian of its rows: the owned degrees of freedom of a cell-centered domain, all
 * degrees of freedom of a box domain, whose partial rows the processes add up. These entries and the entries of
 * forEachRemoteJacobianEntry, summed over all processes and identified by the positions of the degrees of
 * freedom, have to equal the Jacobian on the whole grids up to the error of the numeric differentiation: every
 * entry relative to the largest entry of its row, the entries of the tissue-network blocks relative to the
 * largest entry of these blocks. addRemoteJacobianProduct has to equal the product of these entries with a
 * vector.
 *
 * With Dirichlet conditions on the tissue boundary the problem has a unique solution and is linear in the
 * pressures. With the direct solver MUMPS, Newton's method on the distributed grids has to take as many steps as
 * on the whole grids and reach the same solution, and more steps without the entries of other processes.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#if DUMUX_HAVE_MUMPS
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/mumpssolver.hh>
#include <dumux/multidomain/newtonsolver.hh>
#endif

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
using Position = Dune::FieldVector<double, 3>;
using Communication = Dune::Communication<Dune::MPIHelper::MPICommunicator>;

static_assert(GetPropType<BulkTypeTag, Properties::ModelTraits>::numEq() == 1
              && GetPropType<LowDimTypeTag, Properties::ModelTraits>::numEq() == 1,
              "The test numbers degrees of freedom and assumes one unknown per degree of freedom");

//! The value of a vector at a position, different from the pressure fields
double vectorValue(const Position& x)
{ return 0.7 - 0.4*std::cos(3.0*x[0] - x[2]) + 0.9*x[1]*x[1]; }

//! The coupled problem on two grid views
struct CoupledProblem
{
    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;

    template<class BulkGridView, class LowDimGridView>
    CoupledProblem(const BulkGridView& bulkGridView, const LowDimGridView& lowDimGridView, bool dirichlet)
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
        setValues(sol, [](const auto& x){ return bulkPressure(x, 0); }, [](const auto& x){ return lowDimPressure(x, 0); });

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

    //! Set x to the given functions of the positions of the degrees of freedom
    template<class BulkFunction, class LowDimFunction>
    void setValues(SolutionVector& x, BulkFunction&& bulk, LowDimFunction&& lowDim) const
    {
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
        set(*bulkGridGeometry, x[bulkIdx], bulk);
        set(*lowDimGridGeometry, x[lowDimIdx], lowDim);
    }

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

//! The positions of the degrees of freedom of a grid geometry, by dof index
template<class GridGeometry>
std::vector<Position> dofPositions(const GridGeometry& gridGeometry)
{
    std::vector<Position> positions(gridGeometry.numDofs());
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
            positions[scv.dofIndex()] = scv.dofPosition();
    }
    return positions;
}

//! Global indices of the degrees of freedom: their indices in the whole grids, the network after the tissue
struct GlobalNumbering
{
    explicit GlobalNumbering(const CoupledProblem& whole)
    {
        const std::array<std::vector<Position>, 2> wholePositions{dofPositions(*whole.bulkGridGeometry), dofPositions(*whole.lowDimGridGeometry)};
        std::int64_t offset = 0;
        for (int domain = 0; domain < 2; ++domain)
        {
            for (std::size_t dof = 0; dof < wholePositions[domain].size(); ++dof)
            {
                indices[domain][dofKey(wholePositions[domain][dof])] = offset + static_cast<std::int64_t>(dof);
                positions.push_back(wholePositions[domain][dof]);
            }
            offset += wholePositions[domain].size();
        }
        numBulkDofs = wholePositions[0].size();
    }

    //! the global indices of the degrees of freedom of a coupled problem, by domain and dof index
    std::array<std::vector<std::int64_t>, 2> of(const CoupledProblem& problem) const
    {
        std::array<std::vector<std::int64_t>, 2> result;
        const std::array<std::vector<Position>, 2> localPositions{dofPositions(*problem.bulkGridGeometry), dofPositions(*problem.lowDimGridGeometry)};
        for (int domain = 0; domain < 2; ++domain)
            for (const auto& x : localPositions[domain])
                result[domain].push_back(indices[domain].at(dofKey(x)));
        return result;
    }

    std::array<std::map<DofKey, std::int64_t>, 2> indices;
    std::vector<Position> positions;
    std::int64_t numBulkDofs = 0;
};

//! The rows a process assembles: of interior elements (cell-centered) or of interior and border vertices (box)
template<class GridGeometry>
std::vector<bool> assembledRows(const GridGeometry& gridGeometry)
{
    std::vector<bool> rows(gridGeometry.numDofs(), false);
    const auto& gridView = gridGeometry.gridView();
    if constexpr (GridGeometry::discMethod == DiscretizationMethods::box)
    {
        for (const auto& vertex : vertices(gridView))
            if (vertex.partitionType() == Dune::InteriorEntity || vertex.partitionType() == Dune::BorderEntity)
                rows[gridGeometry.vertexMapper().index(vertex)] = true;
    }
    else
    {
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
            rows[gridGeometry.elementMapper().index(element)] = true;
    }
    return rows;
}

struct Entry
{
    std::int64_t row;
    std::int64_t column;
    double value;
};

using EntryMap = std::map<std::pair<std::int64_t, std::int64_t>, double>;

//! The entries of all processes added up on rank 0
EntryMap sumOnRankZero(const std::vector<Entry>& entries, const Communication& comm)
{
    const int numBytes = static_cast<int>(entries.size()*sizeof(Entry));
    std::vector<int> counts(comm.size()), displacements(comm.size(), 0);
    MPI_Gather(&numBytes, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
    for (int p = 1; p < comm.size(); ++p)
        displacements[p] = displacements[p-1] + counts[p-1];
    std::vector<Entry> all(comm.rank() == 0 ? (displacements.back() + counts.back())/sizeof(Entry) : 0);
    MPI_Gatherv(entries.data(), numBytes, MPI_BYTE, all.data(), counts.data(), displacements.data(), MPI_BYTE, 0, comm);

    EntryMap sum;
    for (const auto& entry : all)
        sum[{entry.row, entry.column}] += entry.value;
    return sum;
}

/*!
 * \brief The entries of the assembled rows of the Jacobian of a coupled problem at its solution
 * \param numRemote the number of entries whose column is a degree of freedom of another process
 * \note Collective.
 */
std::vector<Entry> jacobianEntries(CoupledProblem& problem, const GlobalNumbering& numbering, std::size_t& numRemote)
{
    problem.assembler->assembleJacobianAndResidual(problem.sol);

    const auto indices = numbering.of(problem);
    const std::array<std::vector<bool>, 2> rows{assembledRows(*problem.bulkGridGeometry), assembledRows(*problem.lowDimGridGeometry)};

    std::vector<Entry> entries;
    const auto& jacobian = problem.assembler->jacobian();
    using namespace Dune::Hybrid;
    forEach(std::make_index_sequence<2>{}, [&](auto i)
    {
        forEach(std::make_index_sequence<2>{}, [&](auto j)
        {
            const auto& block = jacobian[i][j];
            for (auto row = block.begin(); row != block.end(); ++row)
                if (rows[i][row.index()])
                    for (auto entry = row->begin(); entry != row->end(); ++entry)
                        entries.push_back({indices[i][row.index()], indices[j][entry.index()], (*entry)[0][0]});
        });
    });

    problem.couplingManager->updateRemoteColumns([&](auto domainIdx, std::size_t dofIdx){ return indices[domainIdx][dofIdx]; });
    numRemote = 0;
    problem.couplingManager->forEachRemoteJacobianEntry([&](auto domainIdx, std::size_t dofIdx, std::size_t eqIdx, std::int64_t column, double value)
    {
        if (rows[domainIdx][dofIdx])
        {
            entries.push_back({indices[domainIdx][dofIdx] + static_cast<std::int64_t>(eqIdx), column, value});
            ++numRemote;
        }
    });

    return entries;
}

/*!
 * \brief The largest difference of an entry relative to a scale, over the entries of both maps
 * \param scale the function scale(row, column) of the sequential entry
 * \param worst a description of the entry with the largest relative difference
 */
template<class Scale>
double maxRelativeDifference(const EntryMap& distributed, const EntryMap& sequential, Scale&& scale, std::string& worst)
{
    double maxRelative = 0.0;
    const auto compare = [&](const auto& key, double d, double s)
    {
        const auto relative = std::abs(d - s)/scale(key.first, key.second);
        if (!std::isfinite(relative) || relative > maxRelative)
        {
            maxRelative = std::isfinite(relative) ? relative : std::numeric_limits<double>::infinity();
            worst = "row " + std::to_string(key.first) + ", column " + std::to_string(key.second)
                    + ": distributed " + std::to_string(d) + ", sequential " + std::to_string(s);
        }
    };
    for (const auto& [key, value] : sequential)
    {
        const auto it = distributed.find(key);
        compare(key, it == distributed.end() ? 0.0 : it->second, value);
    }
    for (const auto& [key, value] : distributed)
        if (!sequential.count(key))
            compare(key, value, 0.0);
    return maxRelative;
}

#if DUMUX_HAVE_MUMPS
template<class LinearSolver>
class StepCountingNewtonSolver : public MultiDomainNewtonSolver<Assembler, LinearSolver, TestCouplingManager>
{
    using ParentType = MultiDomainNewtonSolver<Assembler, LinearSolver, TestCouplingManager>;
public:
    using ParentType::ParentType;

    int numSteps() const
    { return this->numSteps_; }
};

/*!
 * \brief Solve the coupled problem with Newton's method and MUMPS
 * \return the number of Newton steps, -1 if Newton's method did not converge, and the solution
 */
template<class LinearSolver>
std::pair<int, SolutionVector> solve(CoupledProblem& problem, std::shared_ptr<LinearSolver> linearSolver, const Communication& comm)
{
    StepCountingNewtonSolver<LinearSolver> newton(problem.assembler, linearSolver, problem.couplingManager, comm);
    auto x = problem.sol;
    try
    {
        newton.solve(x);
    }
    catch (const NumericalProblem&)
    {
        return {-1, x};
    }
    return {newton.numSteps(), x};
}

//! the direct solver for the distributed grids, optionally with the entries whose column is on another process
auto makeDistributedSolver(const CoupledProblem& problem, bool remoteEntries)
{
    using LinearSolver = DirectSolverMumps<LinearSolverTraits<CoupledProblem::BulkGridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(std::make_tuple(problem.bulkGridGeometry, problem.lowDimGridGeometry));
    if (remoteEntries)
        linearSolver->setAdditionalEntries([couplingManager = problem.couplingManager](bool rebuild, const auto& globalIndex, const auto& add)
        {
            if (rebuild)
                couplingManager->updateRemoteColumns(globalIndex);
            couplingManager->forEachRemoteJacobianEntry([&](auto domainIdx, std::size_t dofIdx, std::size_t eqIdx, std::int64_t column, double value)
            { add(domainIdx, dofIdx, eqIdx, column, value); });
        });
    return linearSolver;
}
#endif

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using BulkGrid = GetPropType<BulkTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimTypeTag, Properties::Grid>;
    const Communication comm(mpiHelper.getCommunicator());
    const Communication self(mpiHelper.getLocalCommunicator());

    int failures = 0;
    const auto fail = [&](const auto& message){ std::cerr << "rank " << comm.rank() << ": " << message << std::endl; ++failures; };

    // the reference on the whole grids with a sequential communicator
    const auto wholeBulkGrid = makeBulkGrid<BulkGrid>(mpiHelper.getLocalCommunicator());
    const auto wholeLowDimGrid = makeNetwork<LowDimGrid>();
    CoupledProblem whole(wholeBulkGrid->leafGridView(), wholeLowDimGrid->leafGridView(), false);
    const GlobalNumbering numbering(whole);

    // the distributed grids
    const auto bulkGrid = makeBulkGrid<BulkGrid>(mpiHelper.getCommunicator());
    auto lowDimGrid = makeNetwork<LowDimGrid>();
    distributeNetwork(*lowDimGrid, bulkGrid->leafGridView());
    CoupledProblem distributed(bulkGrid->leafGridView(), lowDimGrid->leafGridView(), false);

    std::size_t numRemote = 0, unused = 0;
    const auto wholeEntries = jacobianEntries(whole, numbering, unused);
    const auto distributedEntries = jacobianEntries(distributed, numbering, numRemote);
    numRemote = comm.sum(numRemote);

    // the product of the entries whose column is on another process with a vector, at the coupled Jacobian
    {
        auto x = distributed.sol, y = distributed.sol, expected = distributed.sol;
        distributed.setValues(x, [](const auto& p){ return vectorValue(p); }, [](const auto& p){ return vectorValue(p); });
        y = 0.0;
        expected = 0.0;
        distributed.couplingManager->addRemoteJacobianProduct(x, y);
        distributed.couplingManager->forEachRemoteJacobianEntry([&](auto domainIdx, std::size_t dofIdx, std::size_t eqIdx, std::int64_t column, double value)
        { expected[domainIdx][dofIdx][eqIdx] += value*vectorValue(numbering.positions[column]); });

        double difference = 0.0, magnitude = 0.0;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<2>{}, [&](auto i)
        {
            for (std::size_t dof = 0; dof < y[i].size(); ++dof)
            {
                if (!std::isfinite(y[i][dof][0]))
                    difference = std::numeric_limits<double>::infinity();
                difference = std::max(difference, std::abs(y[i][dof][0] - expected[i][dof][0]));
                magnitude = std::max(magnitude, std::abs(expected[i][dof][0]));
            }
        });
        difference = comm.max(difference);
        magnitude = comm.max(magnitude);
        if (comm.rank() == 0)
            std::cout << "product with a vector: max |addRemoteJacobianProduct - entries * x| = " << difference
                      << ", max |entries * x| = " << magnitude << std::endl;
        if (difference > 1e-12*magnitude)
            fail("addRemoteJacobianProduct differs from the product with the entries");
        if (numRemote > 0 && magnitude == 0.0)
            fail("the product with the entries of other processes vanishes");
    }

    const auto wholeJacobian = sumOnRankZero(wholeEntries, self);
    const auto distributedJacobian = sumOnRankZero(distributedEntries, comm);
    if (comm.rank() == 0)
    {
        // numeric differentiation amplifies round-off differences of the residuals, so the entries are compared
        // relative to the largest entry of their row, the tissue-network blocks relative to their largest entry
        std::map<std::int64_t, double> rowMax;
        double couplingMax = 0.0;
        const auto isCoupling = [&](std::int64_t row, std::int64_t column)
        { return (row < numbering.numBulkDofs) != (column < numbering.numBulkDofs); };
        for (const auto& [key, value] : wholeJacobian)
        {
            rowMax[key.first] = std::max(rowMax[key.first], std::abs(value));
            if (isCoupling(key.first, key.second))
                couplingMax = std::max(couplingMax, std::abs(value));
        }

        std::string worstRow, worstCoupling;
        const auto rowDifference = maxRelativeDifference(distributedJacobian, wholeJacobian,
            [&](std::int64_t row, std::int64_t){ return rowMax.at(row); }, worstRow);
        EntryMap wholeCoupling, distributedCoupling;
        for (const auto& [key, value] : wholeJacobian)
            if (isCoupling(key.first, key.second))
                wholeCoupling[key] = value;
        for (const auto& [key, value] : distributedJacobian)
            if (isCoupling(key.first, key.second))
                distributedCoupling[key] = value;
        const auto couplingDifference = maxRelativeDifference(distributedCoupling, wholeCoupling,
            [&](std::int64_t, std::int64_t){ return couplingMax; }, worstCoupling);

        std::cout << "Jacobian: max |distributed - sequential| / max |sequential entry of the row| = " << rowDifference
                  << " (" << worstRow << ")\n"
                  << "tissue-network blocks: max |distributed - sequential| / max |sequential entry| = " << couplingDifference
                  << " (" << worstCoupling << ")\n"
                  << numRemote << " entries with columns of other processes" << std::endl;
        if (couplingMax == 0.0)
            fail("no coupling entries in the sequential Jacobian");
        if (rowDifference > 1e-6)
            fail("the Jacobian differs from the sequential one");
        if (couplingDifference > 1e-6)
            fail("the tissue-network blocks of the Jacobian differ from the sequential ones");
    }
    if (comm.size() > 1 && numRemote == 0)
        fail("no entries with columns of other processes");

#if DUMUX_HAVE_MUMPS
    {
        CoupledProblem wholeDirichlet(wholeBulkGrid->leafGridView(), wholeLowDimGrid->leafGridView(), true);
        using SequentialSolver = DirectSolverMumps<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
        const auto [wholeSteps, wholeSolution] = solve(wholeDirichlet, std::make_shared<SequentialSolver>(), self);

        CoupledProblem distributedDirichlet(bulkGrid->leafGridView(), lowDimGrid->leafGridView(), true);
        const auto [steps, solution] = solve(distributedDirichlet, makeDistributedSolver(distributedDirichlet, true), comm);
        CoupledProblem withoutRemote(bulkGrid->leafGridView(), lowDimGrid->leafGridView(), true);
        const auto [stepsWithoutRemote, unusedSolution] = solve(withoutRemote, makeDistributedSolver(withoutRemote, false), comm);

        const auto wholeIndices = numbering.of(wholeDirichlet);
        const auto indices = numbering.of(distributedDirichlet);
        std::vector<double> reference(numbering.positions.size());
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<2>{}, [&](auto i)
        {
            for (std::size_t dof = 0; dof < wholeSolution[i].size(); ++dof)
                reference[wholeIndices[i][dof]] = wholeSolution[i][dof][0];
        });

        double difference = 0.0, magnitude = 0.0;
        const std::array<std::vector<bool>, 2> rows{assembledRows(*distributedDirichlet.bulkGridGeometry), assembledRows(*distributedDirichlet.lowDimGridGeometry)};
        forEach(std::make_index_sequence<2>{}, [&](auto i)
        {
            for (std::size_t dof = 0; dof < solution[i].size(); ++dof)
                if (rows[i][dof])
                {
                    if (!std::isfinite(solution[i][dof][0]))
                        difference = std::numeric_limits<double>::infinity();
                    difference = std::max(difference, std::abs(solution[i][dof][0] - reference[indices[i][dof]]));
                    magnitude = std::max(magnitude, std::abs(reference[indices[i][dof]]));
                }
        });
        difference = comm.max(difference);
        magnitude = comm.max(magnitude);

        if (comm.rank() == 0)
            std::cout << "Newton steps: " << wholeSteps << " sequential, " << steps << " distributed, "
                      << stepsWithoutRemote << " distributed without the entries of other processes (-1: not converged)\n"
                      << "solution: max |distributed - sequential| / max |sequential| = " << difference/magnitude << std::endl;
        if (wholeSteps < 0 || steps != wholeSteps)
            fail("Newton's method takes a different number of steps on the distributed grids");
        if (difference > 1e-10*magnitude)
            fail("the solution differs from the sequential one");
        if (numRemote > 0 && stepsWithoutRemote >= 0 && stepsWithoutRemote <= steps)
            fail("Newton's method converges as fast without the entries of other processes");
    }
#endif

    failures = comm.sum(failures);
    if (comm.rank() == 0)
        std::cout << (failures == 0 ? "PASSED" : "FAILED") << " on " << comm.size() << " processes" << std::endl;

    return failures == 0 ? 0 : 1;
}
