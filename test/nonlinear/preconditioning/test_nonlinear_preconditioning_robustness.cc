// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \brief Compares the largest time step plain Newton and ASPEN can take on the same problem.
 *
 * Both solvers take a single step of increasing size from the same state, with time-step cutting
 * disabled so that a failure is a failure. The largest size each still converges at is the robustness
 * measure the method is meant to improve.
 *
 * The test asserts only that ASPEN converges wherever Newton does, and reports the rest. Whether
 * nonlinear preconditioning actually buys a larger step is the question being measured, not an
 * assumption to be baked into a pass condition.
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

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TwoPIncompressibleTpfa;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;

    SolutionVector uInit(gridGeometry->numDofs());
    problem->applyInitialSolution(uInit);

    using Index = NonlinearPreconditioning::Partition::Index;
    const auto numElements = gridGeometry->gridView().size(0);
    std::vector<std::vector<Index>> influence(numElements);
    for (Index i = 0; i < static_cast<Index>(numElements); ++i)
        for (const auto& dataJ : gridGeometry->connectivityMap()[i])
            influence[i].push_back(dataJ.globalJ);

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
            auto k = static_cast<std::size_t>(relative*2.0);
            bin = bin*2 + std::min(k, std::size_t{1});
        }
        cores[eIdx] = bin;
    }
    auto partition = std::make_shared<NonlinearPreconditioning::Partition>(influence, cores, 0);

    // Both solvers are held to the same standard: the true global residual must fall below
    // residualTolerance. Comparing Newton's default relative-shift criterion against an absolute
    // residual criterion would not be a like-for-like robustness statement.
    const Scalar residualTolerance = 1e-9;

    const auto globalResidualNorm = [&] (const SolutionVector& u, Scalar dt)
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(u);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uInit);
        assembler->setLinearSystem();
        assembler->assembleResidual(u);
        Scalar norm = 0.0;
        for (const auto& block : assembler->residual())
            norm = std::max(norm, block.infinity_norm());
        return norm;
    };

    const auto runNewton = [&] (Scalar dt, bool lineSearch)
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);

        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uInit);
        auto linearSolver = std::make_shared<LinearSolver>();
        NewtonSolver<Assembler, LinearSolver> newton(assembler, linearSolver);
        newton.setUseLineSearch(lineSearch);
        try {
            if (!newton.apply(u))
                return false;
            return globalResidualNorm(u, dt) < residualTolerance;
        }
        catch (const Dune::Exception&) { return false; }
    };

    const auto runSchwarz = [&] (Scalar dt)
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);

        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);
        NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
            problem, gridGeometry, gridVariables, partition, influence);
        solver.setTimeLoop(timeLoop);
        solver.setPreviousSolution(uInit);
        solver.setLocalResidualTolerance(1e-12);
        solver.setResidualTolerance(residualTolerance);
        solver.setResidualReduction(0.0);   // this comparison is defined by the absolute bar alone
        try {
            const auto report = solver.solve(u);
            if (!report.converged)
                return false;
            if (globalResidualNorm(u, dt) >= residualTolerance)
                return false;

            // a converged residual proves u solves the discrete system, but an unphysical saturation
            // would mean the system was solved on a branch the model does not represent
            for (std::size_t i = 0; i < u.size(); ++i)
                if (u[i][1] < -1e-10 || u[i][1] > 1.0 + 1e-10)
                {
                    std::cerr << "  ASPEN reached an unphysical saturation " << u[i][1]
                              << " at dt = " << dt << std::endl;
                    return false;
                }
            return true;
        }
        catch (const Dune::Exception&) { return false; }
    };

    std::cout << "\n  all three are required to reach a true global residual below "
              << std::scientific << std::setprecision(1) << residualTolerance << "\n\n";
    std::cout << "  time step [s]   Newton      Newton+LS   ASPEN (4 subdomains)\n";
    std::cout << "  ------------------------------------------------------------\n";

    Scalar largestNewton = 0.0, largestLineSearch = 0.0, largestSchwarz = 0.0;
    for (Scalar dt = 250.0; dt <= 4e6; dt *= 2.0)
    {
        const bool newtonOk = runNewton(dt, false);
        const bool lineSearchOk = runNewton(dt, true);
        const bool schwarzOk = runSchwarz(dt);

        if (newtonOk) largestNewton = std::max(largestNewton, dt);
        if (lineSearchOk) largestLineSearch = std::max(largestLineSearch, dt);
        if (schwarzOk) largestSchwarz = std::max(largestSchwarz, dt);

        std::cout << "  " << std::setw(12) << std::scientific << std::setprecision(2) << dt
                  << std::setw(12) << (newtonOk ? "ok" : "FAILED")
                  << std::setw(12) << (lineSearchOk ? "ok" : "FAILED")
                  << std::setw(16) << (schwarzOk ? "ok" : "FAILED") << std::endl;
    }

    std::cout << "\n  largest converging step: Newton " << largestNewton
              << ", Newton+LS " << largestLineSearch
              << ", ASPEN " << largestSchwarz << std::endl;
    const auto bestNewton = std::max(largestNewton, largestLineSearch);
    if (bestNewton > 0.0)
        std::cout << "  ratio ASPEN / best Newton = " << std::fixed << std::setprecision(2)
                  << largestSchwarz/bestNewton << std::endl;

    if (largestSchwarz < bestNewton)
    {
        std::cerr << "FAIL: ASPEN is less robust than plain Newton on this problem" << std::endl;
        return 1;
    }

    std::cout << "\n  ASPEN is at least as robust as plain Newton" << std::endl;
    return 0;
}
