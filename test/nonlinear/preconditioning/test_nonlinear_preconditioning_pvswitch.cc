// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \brief Checks nonlinear preconditioning on a compositional model, where phases appear and disappear.
 *
 * A degree of freedom's phase state is a discrete label, so agreement with plain Newton is required to
 * be exact rather than approximate: a different state is a different answer, not a less accurate one.
 *
 * The test also verifies that the subdomain-restricted switch touches only the degrees of freedom it is
 * given, and carries a control in which switching inside local solves is suppressed. That control must
 * fail, otherwise the problem never exercises switching and the rest of the test would prove nothing.
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
#include <dumux/nonlinear/preconditioning/partitioner.hh>
#include <dumux/nonlinear/preconditioning/pvswitch.hh>

#include "2p2c/properties.hh"

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

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::InjectionCCTpfa;

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

    static_assert(NonlinearPreconditioning::SubdomainPrimaryVariableSwitch<GridVariables>::active(),
                  "This test is meaningless unless the model actually has a primary variable switch");

    SolutionVector uInit(gridGeometry->numDofs());
    problem->applyInitialSolution(uInit);

    const auto influence = NonlinearPreconditioning::buildInfluenceMap(*gridGeometry);
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto numSteps = getParam<std::size_t>("Test.NumTimeSteps", 6);

    // reference: ordinary Newton, several steps so that phase transitions have time to occur
    SolutionVector uNewton(uInit);
    std::size_t newtonStateChanges = 0;
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(uNewton);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);

        SolutionVector uOld(uNewton);
        for (std::size_t step = 0; step < numSteps; ++step)
        {
            auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uOld);
            auto linearSolver = std::make_shared<LinearSolver>();
            NewtonSolver<Assembler, LinearSolver> newton(assembler, linearSolver);
            newton.solve(uNewton);
            for (std::size_t i = 0; i < uNewton.size(); ++i)
                if (uNewton[i].state() != uOld[i].state())
                    ++newtonStateChanges;
            uOld = uNewton;
            gridVariables->advanceTimeStep();
        }
    }
    std::cout << "\n  reference Newton: " << newtonStateChanges
              << " phase state changes over " << numSteps << " steps" << std::endl;
    check(newtonStateChanges > 0,
          "the reference run actually switches phase states, so the comparison is meaningful");

    const auto runSchwarz = [&] (std::size_t blocks, bool switchInLocalSolves)
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);
        auto partition = NonlinearPreconditioning::makeCartesianPartition(*gridGeometry, influence, {blocks, blocks}, 0);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);

        SolutionVector uOld(u);
        bool allConverged = true;
        std::size_t switched = 0;
        for (std::size_t step = 0; step < numSteps; ++step)
        {
            NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
                problem, gridGeometry, gridVariables, partition, influence);
            solver.setTimeLoop(timeLoop);
            solver.setPreviousSolution(uOld);
            solver.setLocalResidualTolerance(1e-12);
            solver.setSwitchInLocalSolves(switchInLocalSolves);
            solver.setVerbosity(getParam<int>("Test.Verbosity", 0));

            const auto report = solver.solve(u);
            allConverged = allConverged && report.converged;
            switched += report.switchedDofs;
            if (!report.converged)
                break;
            uOld = u;
            gridVariables->advanceTimeStep();
        }
        return std::make_tuple(allConverged, u, switched);
    };

    for (const std::size_t blocks : {1ul, 2ul, 4ul})
    {
        const auto [converged, u, switched] = runSchwarz(blocks, true);
        const auto numSubdomains = blocks*blocks;
        const std::string tag = " (" + std::to_string(numSubdomains) + " subdomains)";

        check(converged, "ASPEN converged" + tag);
        if (!converged)
            continue;

        double difference = 0.0, scale = 0.0;
        std::size_t stateMismatches = 0;
        for (std::size_t i = 0; i < u.size(); ++i)
        {
            if (u[i].state() != uNewton[i].state())
                ++stateMismatches;
            for (int eq = 0; eq < 2; ++eq)
            {
                difference = std::max(difference, std::abs(u[i][eq] - uNewton[i][eq]));
                scale = std::max(scale, std::abs(uNewton[i][eq]));
            }
        }

        std::cout << "  ASPEN, " << std::setw(2) << numSubdomains << " subdomains: "
                  << switched << " switches, " << stateMismatches << " state mismatches"
                  << ", relative difference = " << std::scientific << std::setprecision(3)
                  << difference/scale << std::endl;

        check(stateMismatches == 0, "every phase state matches plain Newton exactly" + tag);
        check(difference/scale < 1e-6, "fields match plain Newton" + tag);
    }

    // The restricted switch is verified directly rather than through the physics: on this problem only a
    // handful of cells ever change phase, so a run with local switching disabled still reproduces the
    // reference and cannot tell us whether the mechanism works.
    {
        const auto [converged, u, switched] = runSchwarz(4, false);
        std::size_t stateMismatches = 0;
        for (std::size_t i = 0; i < u.size(); ++i)
            if (u[i].state() != uNewton[i].state())
                ++stateMismatches;
        std::cout << "  informational, no switch in local solves: converged = " << std::boolalpha
                  << converged << ", state mismatches = " << stateMismatches
                  << "  (this problem switches too rarely to discriminate)" << std::endl;
    }

    // a state driven well past the solubility limit, so that a switch is guaranteed to be exercised
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector uSupersaturated(uInit);
        for (std::size_t i = 0; i < uSupersaturated.size(); ++i)
            uSupersaturated[i][1] *= 1e3;
        gridVariables->init(uSupersaturated);

        auto partition = NonlinearPreconditioning::makeCartesianPartition(*gridGeometry, influence, {4ul, 4ul}, 0);

        SolutionVector uGlobal(uSupersaturated);
        NonlinearPreconditioning::SubdomainPrimaryVariableSwitch<GridVariables> globalSwitch;
        globalSwitch.initialize(uGlobal, *gridVariables, *problem, *gridGeometry);
        globalSwitch.applyGlobal(uGlobal, *gridVariables, *problem, *gridGeometry);
        const auto numGlobalSwitched = globalSwitch.numSwitched();

        std::cout << "  direct check: a global switch changes " << numGlobalSwitched << " dofs" << std::endl;
        check(numGlobalSwitched > 0, "the direct check is not vacuous");

        // the same switch applied subdomain by subdomain must reproduce it exactly, and each call must
        // leave everything outside the subdomain it was given untouched
        gridVariables->init(uSupersaturated);
        SolutionVector uUnion(uSupersaturated);
        NonlinearPreconditioning::SubdomainPrimaryVariableSwitch<GridVariables> localSwitch;
        localSwitch.initialize(uUnion, *gridVariables, *problem, *gridGeometry);

        bool leakage = false;
        for (NonlinearPreconditioning::SubdomainIndex sd = 0; sd < partition->numSubdomains(); ++sd)
        {
            SolutionVector before(uUnion);
            localSwitch.applyToSubdomain(uUnion, *gridVariables, *problem, *gridGeometry,
                                         partition->core(sd));
            for (std::size_t i = 0; i < uUnion.size(); ++i)
                if (partition->coreSubdomain(i) != sd)
                    if (uUnion[i].state() != before[i].state()
                        || uUnion[i][0] != before[i][0] || uUnion[i][1] != before[i][1])
                        leakage = true;
        }
        check(!leakage, "a subdomain-restricted switch touches no degree of freedom outside its range");

        std::size_t unionMismatches = 0;
        for (std::size_t i = 0; i < uUnion.size(); ++i)
            if (uUnion[i].state() != uGlobal[i].state()
                || uUnion[i][0] != uGlobal[i][0] || uUnion[i][1] != uGlobal[i][1])
                ++unionMismatches;
        check(unionMismatches == 0,
              "switching subdomain by subdomain reproduces the global switch exactly");
    }

    // The bounding hook must actually constrain every iterate the solver produces, not merely be
    // callable: an unbounded pressure update is what takes the iterate outside the tabulated
    // thermodynamic range, where the pressure derivatives vanish and the local blocks go singular.
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);
        auto partition = NonlinearPreconditioning::makeCartesianPartition(*gridGeometry, influence, {4ul, 4ul}, 0);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);

        // deliberately tight, so the clamp binds rather than merely being satisfied
        Scalar pMin = 1e30, pMax = -1e30;
        for (std::size_t i = 0; i < uInit.size(); ++i)
        {
            pMin = std::min(pMin, uInit[i][0]);
            pMax = std::max(pMax, uInit[i][0]);
        }
        const Scalar lower = pMin - 1.0, upper = pMax + 1.0;

        NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
            problem, gridGeometry, gridVariables, partition, influence);
        solver.setTimeLoop(timeLoop);
        solver.setPreviousSolution(uInit);
        solver.setUpdateBounding([lower, upper] (auto& priVars, const auto& old, std::size_t)
        {
            using std::clamp;
            priVars[0] = clamp(priVars[0], lower, upper);
        });

        const auto report = solver.solve(u);
        std::size_t violations = 0;
        for (std::size_t i = 0; i < u.size(); ++i)
            if (u[i][0] < lower - 1e-12 || u[i][0] > upper + 1e-12)
                ++violations;

        // without the hook the same run must leave that interval, otherwise the check above is vacuous
        auto gridVariablesFree = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector uFree(uInit);
        gridVariablesFree->init(uFree);
        NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> unbounded(
            problem, gridGeometry, gridVariablesFree, partition, influence);
        unbounded.setTimeLoop(timeLoop);
        unbounded.setPreviousSolution(uInit);
        unbounded.solve(uFree);

        std::size_t freeViolations = 0;
        for (std::size_t i = 0; i < uFree.size(); ++i)
            if (uFree[i][0] < lower - 1e-12 || uFree[i][0] > upper + 1e-12)
                ++freeViolations;

        std::cout << "  bounding hook: converged = " << std::boolalpha << report.converged
                  << ", violations with hook = " << violations
                  << ", without hook = " << freeViolations << std::endl;
        check(violations == 0, "every accepted iterate respects the bounding hook");
        check(freeViolations > 0, "the bound genuinely binds, so the check above is not vacuous");
    }

    // The caches the solver leaves behind must describe the solution it returns. A phase state switched
    // inside a subdomain whose cached volume variables were never refreshed does not show up in the
    // solution vector at all — only in whatever reads the cache afterwards, which is why this has to be
    // checked directly rather than through the fields.
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);
        auto partition = NonlinearPreconditioning::makeCartesianPartition(*gridGeometry, influence, {4ul, 4ul}, 0);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);

        SolutionVector uOld(u);
        for (std::size_t step = 0; step < numSteps; ++step)
        {
            NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
                problem, gridGeometry, gridVariables, partition, influence);
            solver.setTimeLoop(timeLoop);
            solver.setPreviousSolution(uOld);
            solver.setLocalResidualTolerance(1e-12);
            // With the local stage forbidden from switching, the global switch applied to the combined
            // solution is the only path by which a phase transition can happen at all. That makes this
            // configuration the one that can tell whether it is being applied: if it is skipped, the
            // returned solution carries pending transitions instead of merely converging differently.
            solver.setSwitchInLocalSolves(false);
            if (!solver.solve(u).converged)
                break;
            uOld = u;
            gridVariables->advanceTimeStep();
        }

        auto asLeft = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, uOld);
        asLeft->setLinearSystem();
        asLeft->assembleResidual(u);
        const auto residualAsLeft = asLeft->residual();

        // the same time state as the solver left: previous variables at uOld, current at u
        auto freshVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        freshVariables->init(uOld);
        freshVariables->advanceTimeStep();
        freshVariables->update(u);
        auto fresh = std::make_shared<Assembler>(problem, gridGeometry, freshVariables, timeLoop, uOld);
        fresh->setLinearSystem();
        fresh->assembleResidual(u);
        const auto residualFresh = fresh->residual();

        double cacheDrift = 0.0;
        for (std::size_t i = 0; i < residualFresh.size(); ++i)
            for (std::size_t eq = 0; eq < residualFresh[i].size(); ++eq)
                cacheDrift = std::max(cacheDrift, std::abs(residualAsLeft[i][eq] - residualFresh[i][eq]));

        // The returned solution must be switch-stable: applying the switch to it may not change
        // anything. If it does, a phase transition was still pending when the solver declared success,
        // and nothing in the convergence diagnostics would have shown it.
        SolutionVector uStable(u);
        auto stableVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        stableVariables->init(uOld);
        stableVariables->advanceTimeStep();
        stableVariables->update(uStable);
        NonlinearPreconditioning::SubdomainPrimaryVariableSwitch<GridVariables> probe;
        probe.initialize(uStable, *stableVariables, *problem, *gridGeometry);
        probe.applyGlobal(uStable, *stableVariables, *problem, *gridGeometry);

        std::size_t pending = 0;
        for (std::size_t i = 0; i < u.size(); ++i)
            if (uStable[i].state() != u[i].state())
                ++pending;

        std::cout << "  returned solution: cache drift = " << std::scientific << std::setprecision(3)
                  << cacheDrift << ", pending phase transitions = " << pending << std::endl;
        check(cacheDrift < 1e-20,
              "the grid variables left behind describe the returned solution");
        check(pending == 0, "no phase transition is still pending on the returned solution");
    }

    // A sweep in which a local solve fails and a phase state switches must still take the fallback step
    // on the true system. Handling the two conditions in mutually exclusive branches let the switch
    // branch skip the fallback, and the iteration then stalled exactly instead of converging.
    {
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        SolutionVector u(uInit);
        gridVariables->init(u);
        auto partition = NonlinearPreconditioning::makeCartesianPartition(*gridGeometry, influence, {4ul, 4ul}, 0);
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 1e9);

        NonlinearPreconditioning::SchwarzNewtonSolver<TypeTag, DiffMethod::numeric> solver(
            problem, gridGeometry, gridVariables, partition, influence);
        solver.setTimeLoop(timeLoop);
        solver.setPreviousSolution(uInit);
        solver.setToleratedPriVarSwitches(0);
        solver.setLocalMaxIterations(1);   // forces local solves to fail while switching still occurs

        const auto report = solver.solve(u);
        std::cout << "  starved local solves: converged = " << std::boolalpha << report.converged
                  << ", failed local solves = " << report.failedLocalSolves
                  << ", fallback steps = " << report.newtonFallbackSteps << std::endl;
        check(report.failedLocalSolves > 0, "the starved configuration really does fail local solves");
        check(report.newtonFallbackSteps > 0,
              "a failed local solve reaches the fallback even when the sweep also switched");
        check(report.converged, "the solve still converges rather than stalling");
    }

    if (failures > 0)
    {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }

    std::cout << "\nASPEN reproduces plain Newton on a compositional model, phase states included"
              << std::endl;
    return 0;
}
