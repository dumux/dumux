// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \brief Wall-clock cost of one nonlinearly preconditioned step against one Newton step.
 *
 * Run with grid variable caching enabled, as production setups are. Reports the cost of a single time
 * step at a size both solvers manage, so that the robustness advantage measured elsewhere can be
 * weighed against what it costs.
 */
#include <config.h>

#include <chrono>
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
#include <dumux/nonlinear/preconditioning/partitioner.hh>

#include "2p/properties.hh"

namespace {

double secondsSince(const std::chrono::steady_clock::time_point& start)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TwoPIncompressibleTpfaCached;

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

    const auto influence = NonlinearPreconditioning::buildInfluenceMap(*gridGeometry);
    const auto dt = getParam<Scalar>("Benchmark.TimeStepSize", 4e3);

    std::cout << "\n  cells: " << gridGeometry->gridView().size(0)
              << ", degrees of freedom: " << gridGeometry->numDofs()
              << ", time step: " << std::scientific << std::setprecision(2) << dt << " s\n\n";

    double newtonTime = 0.0;
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);

        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uInit);
        auto linearSolver = std::make_shared<LinearSolver>();
        NewtonSolver<Assembler, LinearSolver> newton(assembler, linearSolver);

        const auto start = std::chrono::steady_clock::now();
        const bool ok = newton.apply(u);
        newtonTime = secondsSince(start);
        std::cout << "  plain Newton" << std::setw(28) << (ok ? "converged" : "FAILED")
                  << std::setw(12) << std::fixed << std::setprecision(3) << newtonTime << " s" << std::endl;
    }

    for (const std::size_t blocks : {2ul, 4ul, 8ul})
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);

        auto partition = NonlinearPreconditioning::makeCartesianPartition(*gridGeometry, influence, {blocks, blocks}, 0);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);

        NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
            problem, gridGeometry, gridVariables, partition, influence);
        solver.setTimeLoop(timeLoop);
        solver.setPreviousSolution(uInit);
        solver.setLocalResidualTolerance(1e-10);
        solver.setVerbosity(getParam<int>("Benchmark.Verbosity", 0));

        const auto start = std::chrono::steady_clock::now();
        const auto report = solver.solve(u);
        const auto time = secondsSince(start);

        std::cout << "  ASPEN, " << std::setw(3) << partition->numSubdomains() << " subdomains"
                  << std::setw(22) << (report.converged ? "converged" : "FAILED")
                  << std::setw(12) << std::fixed << std::setprecision(3) << time << " s"
                  << "   (" << report.outerIterations << " outer, "
                  << report.globalStages << " global, "
                  << report.localIterations << " local)"
                  << "   relative to Newton: " << std::setprecision(2) << time/newtonTime << "x"
                  << std::endl;
    }

    std::cout << std::endl;
    return 0;
}
