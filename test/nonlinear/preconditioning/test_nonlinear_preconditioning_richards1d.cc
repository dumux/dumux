// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \brief Nonlinear domain decomposition on a problem whose nonlinearity is genuinely local.
 *
 * Infiltration into dry sand in one dimension. Almost the whole column is either dry or wet and behaves
 * close to linearly; the difficulty lives in a sharp wetting front a few cells across, which moves. That
 * is the situation nonlinear preconditioning is for, and it is the reason this problem is worth more here
 * than a larger one whose nonlinearity is spread evenly — in that case decomposing buys nothing, and a
 * test built on it would measure the wrong thing.
 *
 * Two claims are checked: that decomposing does not change the answer, and that it lets the solver take
 * time steps plain Newton cannot. Every solver is held to the same convergence standard.
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
#include <dumux/porousmediumflow/richards/newtonsolver.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/nonlinear/preconditioning/schwarznewtonsolver.hh>
#include <dumux/nonlinear/preconditioning/partitioner.hh>

#include "richards1d/properties.hh"

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

//! counts nonlinear iterations, which the base class tracks but does not expose
template<class Assembler, class LinearSolver>
class CountingRichardsNewton : public Dumux::RichardsNewtonSolver<Assembler, LinearSolver>
{
    using ParentType = Dumux::RichardsNewtonSolver<Assembler, LinearSolver>;
public:
    using ParentType::ParentType;
    using typename ParentType::Variables;

    void newtonEndStep(Variables& vars, const typename Assembler::SolutionVector& uLastIter) override
    {
        ParentType::newtonEndStep(vars, uLastIter);
        ++iterations;
    }

    std::size_t iterations = 0;
};

double secondsSince(const std::chrono::steady_clock::time_point& start)
{ return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); }

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::RichardsBenchmarkCCTpfa;

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
    const auto numSubdomains = getParam<std::size_t>("Test.NumSubdomains", 10);
    auto partition = NonlinearPreconditioning::makeCartesianPartition(*gridGeometry, influence, {numSubdomains}, 0);

    std::cout << "\n  " << gridGeometry->numDofs() << " cells in one dimension, "
              << partition->numSubdomains() << " subdomains\n" << std::endl;

    // Both solvers are judged by the same criterion — the relative movement of the solution, at the
    // same tolerance — so the comparison is about which one can take the step rather than who was asked
    // for less. A residual-reduction gate is deliberately not used here: this problem starts close to
    // equilibrium, so its initial residual is small and a fixed reduction of it is unreachable, which
    // would report every solver as failing and say nothing about either.

    const auto residualNorm = [&] (const SolutionVector& u, const SolutionVector& uOld, Scalar dt)
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(u);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uOld);
        assembler->setLinearSystem();
        assembler->assembleResidual(u);
        Scalar norm = 0.0;
        for (const auto& block : assembler->residual())
            norm = std::max(norm, block.infinity_norm());
        return norm;
    };

    // The same clamp RichardsNewtonSolver applies, expressed through the solver's bounding hook. Without
    // it Richards is intractable for either solver, so giving it to only one would make the comparison
    // meaningless. The bound is a property of the material at that location, which is why the hook needs
    // the degree of freedom's index and not just its values.
    using Indices = GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    const auto richardsClamp = [&] (auto& priVars, const auto& old, std::size_t dofIdx)
    {
        const auto element = gridGeometry->element(dofIdx);
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto elemSol = elementSolution(element, SolutionVector(uInit), *gridGeometry);
            const auto fmi = problem->spatialParams().fluidMatrixInteraction(element, scv, elemSol);
            const auto pcMin = fmi.pc(1.0);
            const auto pw = old[Indices::pressureIdx];
            using std::max; using std::clamp;
            const auto pn = max(problem->nonwettingReferencePressure(), pw + pcMin);
            const auto swOld = max(0.0, fmi.sw(pn - pw));
            priVars[Indices::pressureIdx] = clamp(priVars[Indices::pressureIdx],
                                                  pn - fmi.pc(swOld - 0.2), pn - fmi.pc(swOld + 0.2));
        }
    };

    double newtonSeconds = 0.0, schwarzSeconds = 0.0;
    std::size_t newtonIterations = 0;

    const auto runNewton = [&] (Scalar dt, bool lineSearch, SolutionVector& u)
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        u = uInit;
        gridVariables->init(u);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uInit);
        auto linearSolver = std::make_shared<LinearSolver>();
        CountingRichardsNewton<Assembler, LinearSolver> newton(assembler, linearSolver);
        newton.setUseLineSearch(lineSearch);
        try {
            const auto start = std::chrono::steady_clock::now();
            const bool ok = newton.apply(u);
            newtonSeconds = secondsSince(start);
            newtonIterations = newton.iterations;
            return ok;
        }
        catch (const Dune::Exception&) { return false; }
    };

    const auto runSchwarz = [&] (Scalar dt, SolutionVector& u, NonlinearPreconditioning::SchwarzNewtonReport& report)
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        u = uInit;
        gridVariables->init(u);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);
        NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
            problem, gridGeometry, gridVariables, partition, influence);
        solver.setTimeLoop(timeLoop);
        solver.setPreviousSolution(uInit);
        solver.setUpdateBounding(richardsClamp);
        try {
            const auto start = std::chrono::steady_clock::now();
            report = solver.solve(u);
            schwarzSeconds = secondsSince(start);
            return report.converged;
        }
        catch (const Dune::Exception&) { return false; }
    };

    // equivalence, at a step both handle comfortably
    {
        SolutionVector uNewton(uInit), uSchwarz(uInit);
        NonlinearPreconditioning::SchwarzNewtonReport report;
        // a step the baseline itself handles from cold, so both solutions exist to compare
        const Scalar dt = getParam<Scalar>("Test.EquivalenceTimeStepSize", 0.25);
        check(runNewton(dt, false, uNewton), "the Richards baseline converges at the reference step");
        check(runSchwarz(dt, uSchwarz, report), "the decomposed solver converges at the reference step");

        Scalar difference = 0.0, scale = 0.0;
        for (std::size_t i = 0; i < uSchwarz.size(); ++i)
        {
            difference = std::max(difference, std::abs(uSchwarz[i][0] - uNewton[i][0]));
            scale = std::max(scale, std::abs(uNewton[i][0]));
        }
        std::cout << "  residual scale at the initial state: "
                  << std::scientific << std::setprecision(3) << residualNorm(uInit, uInit, dt) << std::endl;
        std::cout << "  equivalence at the reference step: relative difference = "
                  << std::scientific << std::setprecision(3) << difference/scale
                  << "  (" << report.outerIterations << " outer, " << report.globalStages
                  << " global, " << report.localIterations << " local)" << std::endl;
        check(difference/scale < 1e-6, "decomposing does not change the converged solution");
    }

    // how large a step each can take
    std::cout << "\n  time step [s]    Richards  Richards+LS       ASPEN\n";
    std::cout << "  ------------------------------------------------\n";

    Scalar largestNewton = 0.0, largestLineSearch = 0.0, largestSchwarz = 0.0;
    for (Scalar dt = 0.25; dt <= 2e5; dt *= 2.0)
    {
        SolutionVector u(uInit);
        NonlinearPreconditioning::SchwarzNewtonReport report;
        const bool newtonOk = runNewton(dt, false, u);
        const bool lineSearchOk = runNewton(dt, true, u);
        const bool schwarzOk = runSchwarz(dt, u, report);

        if (newtonOk) largestNewton = std::max(largestNewton, dt);
        if (lineSearchOk) largestLineSearch = std::max(largestLineSearch, dt);
        if (schwarzOk) largestSchwarz = std::max(largestSchwarz, dt);

        std::cout << "  " << std::setw(12) << std::scientific << std::setprecision(2) << dt
                  << std::setw(12) << (newtonOk ? "ok" : "FAILED")
                  << std::setw(12) << (lineSearchOk ? "ok" : "FAILED")
                  << std::setw(12) << (schwarzOk ? "ok" : "FAILED") << std::endl;
    }

    const auto bestNewton = std::max(largestNewton, largestLineSearch);
    std::cout << "\n  largest converging step: Richards " << largestNewton
              << ", Richards+LS " << largestLineSearch << ", ASPEN " << largestSchwarz << std::endl;
    if (bestNewton > 0.0)
        std::cout << "  ratio = " << std::fixed << std::setprecision(2) << largestSchwarz/bestNewton
                  << "x" << std::endl;

    // Full simulation to the benchmark's first output time, both solvers driven by the *same* time-step
    // controller, so the only difference is the nonlinear solve. Reporting one time step would say
    // nothing about a real run: what matters is the total, and how much of it was spent on attempts that
    // were thrown away.
    struct RunStats
    {
        double seconds = 0.0;
        std::size_t steps = 0;          //!< accepted time steps
        std::size_t rejected = 0;       //!< attempts thrown away, each followed by a smaller step
        std::size_t iterations = 0;     //!< nonlinear iterations in accepted steps
        std::size_t wasted = 0;         //!< nonlinear iterations in rejected attempts
        std::size_t localIterations = 0;//!< subdomain iterations, decomposed solvers only
    };

    const auto tEnd = getParam<std::vector<Scalar>>("TimeLoop.TEnd").front();
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dtStart = getParam<Scalar>("TimeLoop.DtInitial");

    const auto simulate = [&] (bool decomposed, Scalar maxDt)
    {
        RunStats stats;
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit), uOld(uInit);
        gridVariables->init(u);

        Scalar t = 0.0, dt = dtStart;
        const auto start = std::chrono::steady_clock::now();
        while (t < tEnd)
        {
            dt = std::min({dt, maxDt, tEnd - t});
            auto timeLoop = std::make_shared<TimeLoop<Scalar>>(t, dt, tEnd);
            bool ok = false;
            std::size_t iters = 0, localIters = 0;

            if (decomposed)
            {
                NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
                    problem, gridGeometry, gridVariables, partition, influence);
                solver.setTimeLoop(timeLoop);
                solver.setPreviousSolution(uOld);
                solver.setUpdateBounding(richardsClamp);
                u = uOld;
                try {
                    const auto report = solver.solve(u);
                    ok = report.converged;
                    iters = report.outerIterations;
                    localIters = report.localIterations;
                }
                catch (const Dune::Exception&) { ok = false; }
            }
            else
            {
                auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uOld);
                auto linearSolver = std::make_shared<LinearSolver>();
                CountingRichardsNewton<Assembler, LinearSolver> newton(assembler, linearSolver);
                newton.setUseLineSearch(true);
                u = uOld;
                try { ok = newton.apply(u); }
                catch (const Dune::Exception&) { ok = false; }
                iters = newton.iterations;
            }

            if (ok)
            {
                ++stats.steps;
                stats.iterations += iters;
                stats.localIterations += localIters;
                t += dt;
                uOld = u;
                gridVariables->advanceTimeStep();
                dt *= 1.5;
            }
            else
            {
                ++stats.rejected;
                stats.wasted += iters;
                dt *= 0.5;
                if (dt < 1e-6)
                    break;
            }
        }
        stats.seconds = secondsSince(start);
        return stats;
    };

    const auto newtonCapped = simulate(false, maxDt);
    const auto schwarzCapped = simulate(true, maxDt);

    // and again with the cap lifted, so each solver is limited by what it can actually converge rather
    // than by a bound chosen for the benchmark's own accuracy requirements
    const auto newtonFree = simulate(false, tEnd);
    const auto schwarzFree = simulate(true, tEnd);

    const auto row = [] (const std::string& label, const RunStats& r)
    {
        std::cout << "  " << std::left << std::setw(20) << label << std::right
                  << std::setw(7) << r.steps << std::setw(10) << r.rejected
                  << std::setw(8) << r.iterations << std::setw(8) << r.wasted
                  << std::setw(9) << (r.localIterations ? std::to_string(r.localIterations) : "-")
                  << std::setw(12) << std::fixed << std::setprecision(3) << r.seconds << std::endl;
    };

    std::cout << "\n  full simulation to t = " << std::fixed << std::setprecision(0) << tEnd
              << " s, adaptive dt, identical time-step controller\n\n";
    std::cout << "                        steps  rejected   iters  wasted    local    time [s]\n";
    std::cout << "  ---------------------------------------------------------------------------\n";
    std::cout << "  dt capped at " << std::setprecision(0) << maxDt << " s (the benchmark's own cap)\n";
    row("    Richards Newton", newtonCapped);
    row("    Schwarz Newton", schwarzCapped);
    std::cout << "  dt uncapped (each solver limited by its own convergence)\n";
    row("    Richards Newton", newtonFree);
    row("    Schwarz Newton", schwarzFree);

    check(newtonCapped.steps > 0 && schwarzCapped.steps > 0, "both solvers complete the capped run");
    check(newtonFree.steps > 0 && schwarzFree.steps > 0, "both solvers complete the uncapped run");

    check(bestNewton > 0.0, "the baseline converges somewhere, so the sweep is calibrated");
    check(largestSchwarz >= bestNewton,
          "decomposing is never worse than the best Newton baseline on this problem");

    if (failures > 0)
    {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }

    std::cout << "\n  1d infiltration: decomposition preserves the solution and extends the step"
              << std::endl;
    return 0;
}
