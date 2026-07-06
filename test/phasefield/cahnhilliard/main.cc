// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CahnHilliardTests
 * \brief Test for the Cahn-Hilliard model discretized with PQ2 (hybrid CVFE) elements
 *        for both the concentration and the chemical potential.
 */
#include <config.h>

#include <cmath>
#include <vector>
#include <random>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/discretization/cvfe/interpolate.hh>
#include <dumux/phasefield/common/integrate.hh>
#include <dumux/phasefield/cahnhilliard/flux.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/io/chrono.hh>
#include <dumux/io/format.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#if DUMUX_HAVE_GRIDFORMAT
#include <dumux/io/gridwriter.hh>
#endif

#include "properties.hh"

struct MinScatter
{
    template<class A, class B>
    static void apply(A& a, const B& b)
    { a[0] = std::min(a[0], b[0]); }
};

// create the random initial solution
template<class SolutionVector, class GridGeometry>
SolutionVector createInitialSolution(const GridGeometry& gg)
{
    SolutionVector sol(gg.numDofs());

    // Generate random number and add processor offset
    // For sequential runs `rank` always returns `0`.
    std::mt19937 gen(0); // seed is 0 for deterministic results
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int n = 0; n < sol.size(); ++n)
    {
        sol[n][0] = 0.42 + 0.02*(0.5-dis(gen)) + gg.gridView().comm().rank();
        sol[n][1] = 0.0;
    }

    // We take the value of the processor with the minimum rank
    // and subtract the rank offset
    if (gg.gridView().comm().size() > 1)
    {
        std::bitset<GridGeometry::GridView::dimension+1> activeCodims;
        activeCodims.set(GridGeometry::GridView::dimension); activeCodims.set(GridGeometry::GridView::dimension-1);
        Dumux::MultiCodimVectorCommDataHandle<
            typename GridGeometry::DofMapper,
            SolutionVector,
            GridGeometry::GridView::dimension,
            MinScatter
        > minHandle(gg.dofMapper(), sol, activeCodims);
        gg.gridView().communicate(minHandle, Dune::All_All_Interface, Dune::ForwardCommunication);

        // Remove processor offset
        for (int n = 0; n < sol.size(); ++n)
            sol[n][0] -= std::floor(sol[n][0]);
    }
    return sol;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // We initialize MPI and/or multithreading backend.
    Dumux::initialize(argc, argv);

    // parse the parameter tree (params.input)
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::TestCahnHilliardPQ2;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    // create the grid, grid geometry, and problem
    GridManager<Grid> gridManager;
    gridManager.init();

    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);

    // create the initial solution: a small, deterministic (position-based) perturbation
    // around a mean concentration; the chemical potential starts out at zero everywhere.
    // Being a pure function of position (rather than a stateful RNG advanced dof-by-dof),
    // this is trivially consistent under any grid partitioning.
    SolutionVector sol = createInitialSolution<SolutionVector>(*gridGeometry);

    // the convex-concave (Eyre) split reaction term in problem->source() needs
    // the previous time step's solution; seed it with the initial condition
    problem->updateOldSolution(sol);

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // Under the homogeneous Neumann (no-flux) boundary conditions used by this
    // test, the total concentration is a conserved quantity of the continuous
    // Cahn-Hilliard equation. The quantity the discrete scheme actually
    // conserves exactly (up to Newton tolerance) is the control-volume
    // (vertex) dofs' own lumped mass: their sub-control-volume-face flux
    // forms a closed, self-telescoping subsystem independent of the
    // non-CV (edge/center) dofs, so we deliberately check only that (see
    // PhaseField::integrateControlVolumeField), not the full hybrid field.
    auto totalConcentration = [&](const SolutionVector& s)
    {
        return PhaseField::integrateControlVolumeField(*gridGeometry, gridVariables->curGridVars(), s,
            [](const auto& vars) { return vars.concentration(); });
    };
    const auto initialMass = totalConcentration(sol);

    // The Cahn-Hilliard equation is a gradient flow of the free energy
    // F[c] = integral[f(c) + (gamma/2)|grad c|^2], so F must be
    // non-increasing along the discrete trajectory; a small tolerance
    // accommodates roundoff/quadrature noise around a true decrease.
    auto totalEnergy = [&](const SolutionVector& s)
    {
        return PhaseField::integrateField(*gridGeometry, gridVariables->curGridVars(), s,
            [&](const auto& fvGeometry, const auto& elemVars, const auto& ipCache)
            {
                const CahnHilliardFluxFunctionContext context(fvGeometry, elemVars, ipCache);
                const auto& gradC = context.gradConcentration();
                return problem->freeEnergyDensity(context.concentration()) + 0.5*problem->surfaceTension()*(gradC*gradC);
            });
    };
    const auto initialEnergy = totalEnergy(sol);
    auto previousEnergy = initialEnergy;
    static constexpr Scalar energyIncreaseTolerance = 1e-8;

    // instantiate the time loop using start/end time and initial time step size from params.input
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        Chrono::toSeconds(getParam("TimeLoop.TStart", "0")),
        Chrono::toSeconds(getParam("TimeLoop.InitialTimeStepSize")),
        Chrono::toSeconds(getParam("TimeLoop.TEnd"))
    );
    timeLoop->setMaxTimeStepSize(Chrono::toSeconds(getParam("TimeLoop.MaxTimeStepSize")));

    // assembler, linear solver, and Newton solver
    using Assembler = Dumux::Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = SSORBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    auto oldSol = sol;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    // write the concentration and chemical potential fields, at the actual
    // (order-2) Lagrange resolution of the PQ2 discretization, using GridFormat.
    // Output.Every controls how many timesteps pass between writes (1 = every step).
#if DUMUX_HAVE_GRIDFORMAT
    const int outputEvery = getParam<int>("Output.Every", 1);
    std::vector<Scalar> concentration(sol.size());
    std::vector<Scalar> chemicalPotential(sol.size());
    IO::GridWriter<typename GridGeometry::GridView, 2> vtkWriter{IO::Format::vtu, gridGeometry->gridView(), IO::order<2>};
    auto writeVtk = [&]()
    {
        for (std::size_t i = 0; i < sol.size(); ++i)
        {
            concentration[i] = sol[i][Indices::concentrationIdx];
            chemicalPotential[i] = sol[i][Indices::chemicalPotentialIdx];
        }
        vtkWriter.setPointField("c", concentration);
        vtkWriter.setPointField("mu", chemicalPotential);
        vtkWriter.write(problem->name() + "-" + Fmt::format("{:05d}", timeLoop->timeStepIndex()));
    };
    writeVtk();
#endif

    // the time loop
    timeLoop->start(); do
    {
        solver.solve(sol, *timeLoop);

        const auto currentEnergy = totalEnergy(sol);
        if (currentEnergy > previousEnergy + energyIncreaseTolerance*std::abs(previousEnergy))
            DUNE_THROW(Dune::Exception, "Cahn-Hilliard free energy increased: " << previousEnergy
                                         << " -> " << currentEnergy << " at time " << timeLoop->time());
        previousEnergy = currentEnergy;

        oldSol = sol;
        problem->updateOldSolution(sol);
        gridVariables->advanceTimeStep();

        timeLoop->advanceTimeStep();
        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(solver.suggestTimeStepSize(timeLoop->timeStepSize()));

#if DUMUX_HAVE_GRIDFORMAT
        if (timeLoop->timeStepIndex() % outputEvery == 0 || timeLoop->finished())
            writeVtk();
#endif

    } while (!timeLoop->finished());
    timeLoop->finalize(gridGeometry->gridView().comm());

    const auto finalMass = totalConcentration(sol);
    const auto relativeMassDrift = std::abs(finalMass - initialMass)/initialMass;
    std::cout << "Mass conservation check: initial = " << initialMass << ", final = " << finalMass
              << ", relative drift = " << relativeMassDrift << std::endl;

    // The control-volume dofs' lumped mass is conserved by construction
    // (their scvf fluxes telescope to exactly the boundary flux, independent
    // of grid resolution), so the only remaining drift is Newton's own
    // nonlinear solver tolerance -- observed at ~1e-14/~1e-15 on both the
    // default (100x100) and coarser CI (10x10) grids.
    static constexpr Scalar relativeMassTolerance = 1e-8;
    if (relativeMassDrift > relativeMassTolerance)
        DUNE_THROW(Dune::Exception, "Cahn-Hilliard mass conservation violated: relative drift = "
                                     << relativeMassDrift << " > tolerance " << relativeMassTolerance);

    std::cout << "Energy dissipation check passed: free energy decreased monotonically from "
              << initialEnergy << " to " << previousEnergy << " (within tolerance)" << std::endl;

    return 0;
}
