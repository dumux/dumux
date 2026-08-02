// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \brief Checks that nonlinear preconditioning converges to the same solution as plain Newton.
 *
 * A single time step of the two-phase problem is solved twice, once with an ordinary Newton solver and
 * once with the ASPEN solver, and the resulting fields are compared. The degenerate single-subdomain
 * case is checked separately: there the local stage already solves the whole problem, so the global
 * stage must never run.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/nonlinear/preconditioning/schwarznewtonsolver.hh>
#include <dumux/nonlinear/preconditioning/partition.hh>

#include "2p/properties.hh"

namespace {

int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << what << std::endl;
        ++failures;
    }
}

} // end anonymous namespace

template<class TypeTag>
void runFor(const std::string& label)
{
    using namespace Dumux;
    std::cout << "\n  --- " << label << " ---" << std::endl;

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector uInit(gridGeometry->numDofs());
    problem->applyInitialSolution(uInit);

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, 250.0, 3000.0);

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // reference: an ordinary Newton solve of one time step
    SolutionVector uNewton(uInit);
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(uNewton);

        using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uInit);
        using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
        auto linearSolver = std::make_shared<LinearSolver>();
        NewtonSolver<Assembler, LinearSolver> newton(assembler, linearSolver);
        newton.solve(uNewton);
    }

    using Index = NonlinearPreconditioning::Partition::Index;
    const auto numElements = gridGeometry->gridView().size(0);

    std::vector<std::vector<Index>> influence(numElements);
    for (Index i = 0; i < static_cast<Index>(numElements); ++i)
        for (const auto& dataJ : gridGeometry->connectivityMap()[i])
            influence[i].push_back(dataJ.globalJ);

    const auto runSchwarz = [&] (std::size_t numSubdomainsPerDirection)
    {
        std::vector<NonlinearPreconditioning::SubdomainIndex> cores(numElements);
        const auto& bBoxMin = gridGeometry->bBoxMin();
        const auto& bBoxMax = gridGeometry->bBoxMax();
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            const auto center = element.geometry().center();
            std::size_t bin = 0;
            for (int dir = 0; dir < 2; ++dir)
            {
                const auto relative = (center[dir] - bBoxMin[dir])/(bBoxMax[dir] - bBoxMin[dir]);
                auto k = static_cast<std::size_t>(relative*double(numSubdomainsPerDirection));
                k = std::min(k, numSubdomainsPerDirection - 1);
                bin = bin*numSubdomainsPerDirection + k;
            }
            cores[eIdx] = bin;
        }

        auto partition = std::make_shared<NonlinearPreconditioning::Partition>(influence, cores, 0);

        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);

        NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
            problem, gridGeometry, gridVariables, partition, influence);
        solver.setTimeLoop(timeLoop);
        solver.setPreviousSolution(uInit);
        solver.setLocalResidualTolerance(1e-12);

        const auto report = solver.solve(u);
        return std::make_pair(report, u);
    };

    std::cout << std::scientific << std::setprecision(3);
    for (const std::size_t perDirection : {1ul, 2ul, 4ul})
    {
        const auto [report, u] = runSchwarz(perDirection);
        const auto numSubdomains = perDirection*perDirection;
        const std::string tag = " (" + std::to_string(numSubdomains) + " subdomains)";

        check(report.converged, "ASPEN converged" + tag);

        double difference = 0.0, scale = 0.0;
        for (std::size_t i = 0; i < u.size(); ++i)
            for (int eq = 0; eq < 2; ++eq)
            {
                difference = std::max(difference, std::abs(u[i][eq] - uNewton[i][eq]));
                scale = std::max(scale, std::abs(uNewton[i][eq]));
            }

        std::cout << "  ASPEN, " << numSubdomains << " subdomains: "
                  << report.outerIterations << " outer, "
                  << report.globalStages << " global, "
                  << report.localIterations << " local iterations"
                  << ", relative difference to Newton = " << difference/scale << std::endl;

        check(difference/scale < 1e-6, "ASPEN reproduces the Newton solution" + tag);

        // with a single subdomain the local stage solves the whole problem, so no global solve is needed
        if (numSubdomains == 1)
            check(report.globalStages == 0, "the single-subdomain case never enters the global stage");
        else
            check(report.globalStages > 0, "the multi-subdomain case exercises the global stage" + tag);
    }

}

/*!
 * A convergence gate expressed as an absolute residual threshold is only meaningful once the residual's
 * scale is known, and that scale is a property of the problem: cell size, source magnitude, time step.
 * Where it sits below the threshold, such a gate accepts on the first sweep with the solution untouched
 * and reports success. This checks that the solver measures its own scale instead, by driving the
 * residual scale far down and requiring the solve to still do real work.
 */
template<class TypeTag>
void checkScaleFreeConvergence()
{
    using namespace Dumux;
    std::cout << "\n  --- scale-free convergence ---" << std::endl;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    SolutionVector uInit(gridGeometry->numDofs());
    problem->applyInitialSolution(uInit);

    // a very small step makes the residual of the initial guess tiny in absolute terms, the regime in
    // which a fixed absolute threshold silently becomes no threshold at all
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, 1e-3, 1e9);
    const auto influence = NonlinearPreconditioning::buildInfluenceMap(*gridGeometry);
    auto partition = NonlinearPreconditioning::makeCartesianPartition(*gridGeometry, influence, {2ul, 2ul}, 0);

    const auto runWithReduction = [&] (double reduction)
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);

        NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
            problem, gridGeometry, gridVariables, partition, influence);
        solver.setTimeLoop(timeLoop);
        solver.setPreviousSolution(uInit);
        solver.setResidualReduction(reduction);
        return solver.solve(u);
    };

    const auto loose = runWithReduction(1e-2);
    const auto strict = runWithReduction(1e-12);

    std::cout << "  reduction 1e-2  : converged = " << std::boolalpha << loose.converged
              << ", residual = " << std::scientific << std::setprecision(3) << loose.residualNorm
              << ", outer = " << loose.outerIterations << std::endl;
    std::cout << "  reduction 1e-12 : converged = " << std::boolalpha << strict.converged
              << ", residual = " << std::scientific << std::setprecision(3) << strict.residualNorm
              << ", outer = " << strict.outerIterations << std::endl;

    check(loose.converged && strict.converged, "both reductions converge");

    // If the gate were a fixed absolute threshold it would stop at the same place regardless of what
    // reduction was asked for, which is precisely what makes such a gate silently ineffective on a
    // problem whose residual scale happens to sit below it.
    check(strict.residualNorm < loose.residualNorm,
          "the accepted residual responds to the requested reduction");

    // A reduction target below the floating-point floor of this problem's residual is unreachable by
    // construction. The solve must still terminate, via the settled-solution criterion, rather than
    // spending its whole budget failing to reach it — that is what a fine grid with a small source
    // term produces, and what made the default unusable on a real problem.
    const auto unreachable = runWithReduction(1e-30);
    std::cout << "  reduction 1e-30 : converged = " << std::boolalpha << unreachable.converged
              << ", residual = " << std::scientific << std::setprecision(3) << unreachable.residualNorm
              << ", outer = " << unreachable.outerIterations << std::endl;
    check(unreachable.converged,
          "an unreachable reduction target still terminates, via the settled solution");
    check(unreachable.outerIterations < 30,
          "and does so without exhausting the iteration budget");

    // Note on what this does and does not establish: it shows the reduction criterion is wired and
    // governs where the solve stops. It cannot show that the *default* is scale-free, because on this
    // problem the residual scale is far above any fixed threshold and the outer iteration converges
    // quadratically past whatever bar is set, landing near machine precision either way. A problem whose
    // residuals sit below the threshold is what separates the two, and this is not one.

}

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);
    Parameters::init(argc, argv);

    runFor<Properties::TTag::TwoPIncompressibleTpfa>("grid variable caching disabled");
    runFor<Properties::TTag::TwoPIncompressibleTpfaCached>("grid variable caching enabled");
    checkScaleFreeConvergence<Properties::TTag::TwoPIncompressibleTpfaCached>();

    if (failures > 0)
    {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }

    std::cout << "\nASPEN reproduces the Newton solution on every partition" << std::endl;
    return 0;
}
