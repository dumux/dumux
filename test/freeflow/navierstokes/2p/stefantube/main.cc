// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Stefan-tube slit-capillary evaporation benchmark.
 */

#include <config.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>

#include "properties.hh"
#include "../adapt.hh"

namespace {

template<class Scalar>
struct StefanTubeQoI
{
    Scalar gasLength = 0.0;
    Scalar analyticalGasLength = 0.0;
    Scalar gasLengthError = 0.0;
    Scalar liquidMass = 0.0;
    Scalar evaporationRate = 0.0;
    Scalar analyticalEvaporationRate = 0.0;
    Scalar vaporL2Error = 0.0;
};

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomentumTypeTag = Properties::TTag::TYPETAG_MOMENTUM;
    using MassTypeTag = Properties::TTag::TYPETAG_MASS;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);
    Dune::Timer timer;

    using Grid = GetPropType<MassTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using Scalar = typename Traits::Scalar;
    using SolutionVector = typename Traits::SolutionVector;
    using MassIndices = typename GetPropType<MassTypeTag, Properties::ModelTraits>::Indices;

    const auto tStart = getParam<Scalar>("TimeLoop.TStart", 0.0);
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(tStart, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;

    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());

    std::cout << "Total number of dofs: "
              << massGridGeometry->numDofs()
                   + momentumGridGeometry->numDofs()*MomentumGridGeometry::GridView::dimension
              << std::endl;

    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    auto xOld = x;

    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    const Scalar refineTol = getParam<Scalar>("Adaptive.RefineTolerance");
    const Scalar coarsenTol = getParam<Scalar>("Adaptive.CoarsenTolerance");
    TwoPhaseCahnHilliardGridAdaptIndicator<MassTypeTag> indicator(massGridGeometry);

    using GridDataTransfer = TwoDomainOneGridDataTransfer<
        Grid,
        TwoPhaseCahnHilliardMassGridDataTransfer<MassTypeTag>,
        TwoPhaseCahnHilliardVelocityGridDataTransfer<MomentumTypeTag>
    >;
    GridDataTransfer dataTransfer(
        std::make_shared<TwoPhaseCahnHilliardMassGridDataTransfer<MassTypeTag>>(
            massProblem, massGridGeometry, massGridVariables, x[massIdx]),
        std::make_shared<TwoPhaseCahnHilliardVelocityGridDataTransfer<MomentumTypeTag>>(
            momentumProblem, momentumGridGeometry, momentumGridVariables, x[momentumIdx])
    );

    // Declared here (rather than just before the real time loop) so the optional
    // corner-relaxation diagnostic below can reuse the same assembler/solver types.
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;

    // Refine (and re-sample the analytical initial condition on the new grid, see
    // below) until the indicator marks nothing more to refine -- i.e. until the
    // initial mesh has genuinely converged, subject only to Adaptive.MaxLevel/
    // MinElementSize. A fixed pass count (the previous Adaptive.InitMaxLevel) can
    // stop short of that fixed point if convergence needs more passes than
    // guessed, silently under-resolving the initial mesh; markElements() already
    // returns false (breaking the loop) once nothing changes, so the only reason
    // to bound the loop at all is as a safety net against a runaway indicator.
    static constexpr std::size_t maxInitPasses = 100;
    for (std::size_t i = 0; i < maxInitPasses; ++i)
    {
        indicator.calculate(x[massIdx], refineTol, coarsenTol);
        if (!markElements(gridManager.grid(), indicator))
            break;
        if (!adapt(gridManager.grid(), dataTransfer))
            break;

        xOld = x;
        couplingManager->updateSolution(x);
        massGridVariables->updateAfterGridAdaption(x[massIdx]);
        momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);

        // Initial AMR decisions must always be based on the analytical profile
        // on the current grid. Otherwise transfer artifacts from the previous
        // grid can trigger refinement far away from the diffuse interface.
        momentumProblem->applyInitialSolution(x[momentumIdx]);
        massProblem->applyInitialSolution(x[massIdx]);
        xOld = x;
        couplingManager->updateSolution(x);
        massGridVariables->updateAfterGridAdaption(x[massIdx]);
        momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);
    }

    // Re-sample once more on the final initial grid. This is cheap and keeps
    // the behavior identical when the initial AMR loop exits without adapting.
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    xOld = x;

    couplingManager->init(momentumProblem, massProblem,
                          std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->updateAfterGridAdaption(x[massIdx]);
    momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);

    // --- Corner-wicking diagnostic (off by default) -------------------------
    // Concus & Finn (1969, PNAS 63(2):292-299) prove that for a wedge of
    // half-angle alpha and contact angle theta, NO bounded equilibrium
    // meniscus exists when theta + alpha < 90 deg -- here, theta < 45 deg for
    // our square tube corner (alpha=45 deg). At the shipped ContactAngle
    // (30 deg, 2D and 3D alike) this benchmark is BELOW that threshold: the
    // true physical meniscus wicks the corner without bound (liquid volume
    // ~ O(r), surface height ~ O(1/r) approaching the vertex) -- there is no
    // finite shape to relax toward, so a relaxation phase here is NOT
    // expected to converge. This block exists only to OBSERVE that trend
    // (does the liquid-like region keep advancing into the corner rather
    // than settling, confirming the theorem numerically) -- it does not seed
    // the real run unless Problem.UseRelaxedCornerIC is explicitly set, since
    // whatever it stops at is a mesh/dt-dependent artifact of where the
    // safety cap or discretization happened to catch it, not a physical answer.
    if (getParam<bool>("Problem.EnableCornerRelaxation", false))
    {
        const auto xBeforeRelax = x;
        massProblem->setRelaxationPhase(true);

        auto relaxTimeLoop = std::make_shared<TimeLoop<Scalar>>(Scalar(0.0), dt, Scalar(1e6));
        relaxTimeLoop->setMaxTimeStepSize(maxDt);
        auto relaxAssembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                           std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                           std::make_tuple(momentumGridVariables, massGridVariables),
                                                           couplingManager, relaxTimeLoop, xOld);
        NewtonSolver relaxSolver(relaxAssembler, linearSolver, couplingManager);

        const auto& bMax = massGridGeometry->bBoxMax();
        const Scalar cornerTol = 2.0*getParam<Scalar>("Adaptive.MinElementSize", 0.01);
        auto cornerProtrusion = [&]()
        {
            Scalar maxLiquidX = std::numeric_limits<Scalar>::lowest();
            auto fvGeometry = localView(*massGridGeometry);
            for (const auto& element : elements(massGridGeometry->gridView()))
            {
                fvGeometry.bind(element);
                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto& pos = scv.dofPosition();
                    bool nearCorner = true;
                    for (int dir = 1; dir < MassGridGeometry::GridView::dimensionworld; ++dir)
                        nearCorner = nearCorner && (pos[dir] > bMax[dir] - cornerTol);
                    if (nearCorner && x[massIdx][scv.dofIndex()][MassIndices::phaseFieldIdx] > 0.0)
                        maxLiquidX = std::max(maxLiquidX, pos[0]);
                }
            }
            return maxLiquidX;
        };

        const Scalar relaxConvTol = getParam<Scalar>("Problem.RelaxationConvergenceTolerance", 1e-5);
        const auto maxRelaxSteps = getParam<std::size_t>("Problem.MaxRelaxationSteps", 200);
        std::cout << "\033[1;33m[relax] corner-wicking diagnostic: initial protrusion x="
                  << cornerProtrusion() << "\033[0m" << std::endl;

        relaxTimeLoop->start();
        std::size_t step = 0;
        for (; step < maxRelaxSteps; ++step)
        {
            relaxSolver.solve(x, *relaxTimeLoop);

            Scalar maxDeltaPhi = 0.0;
            for (std::size_t i = 0; i < x[massIdx].size(); ++i)
                maxDeltaPhi = std::max(maxDeltaPhi, std::abs(
                    x[massIdx][i][MassIndices::phaseFieldIdx] - xOld[massIdx][i][MassIndices::phaseFieldIdx]));

            xOld = x;
            momentumGridVariables->advanceTimeStep();
            massGridVariables->advanceTimeStep();
            relaxTimeLoop->advanceTimeStep();
            relaxTimeLoop->setTimeStepSize(relaxSolver.suggestTimeStepSize(relaxTimeLoop->timeStepSize()));

            std::cout << "\033[1;33m[relax] step " << step << " max|dphi|=" << maxDeltaPhi
                      << " corner protrusion x=" << cornerProtrusion() << "\033[0m" << std::endl;

            if (maxDeltaPhi < relaxConvTol)
            {
                std::cout << "\033[1;33m[relax] converged after " << (step + 1) << " steps\033[0m" << std::endl;
                break;
            }
        }
        if (step == maxRelaxSteps)
            std::cout << "\033[1;33m[relax] hit the safety cap of " << maxRelaxSteps
                      << " steps without converging -- consistent with the Concus-Finn"
                         " unbounded-corner prediction at this ContactAngle\033[0m" << std::endl;

        massProblem->setRelaxationPhase(false);

        if (!getParam<bool>("Problem.UseRelaxedCornerIC", false))
        {
            x = xBeforeRelax;
            xOld = xBeforeRelax;
            couplingManager->updateSolution(x);
            massGridVariables->updateAfterGridAdaption(x[massIdx]);
            momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);
            std::cout << "\033[1;33m[relax] diagnostic only -- reverted to the original IC for the"
                         " real run (set Problem.UseRelaxedCornerIC=true to keep the relaxed state)"
                         "\033[0m" << std::endl;
        }
    }

    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    auto computeQoI = [&](const Scalar time)
    {
        StefanTubeQoI<Scalar> qoi;

        Scalar liquidVolume = 0.0;
        Scalar evaporationMass = 0.0;
        Scalar vaporErr2 = 0.0;
        Scalar vaporNorm2 = 0.0;

        auto fvGeometry = localView(*massGridGeometry);
        for (const auto& element : elements(massGridGeometry->gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& priVars = x[massIdx][scv.dofIndex()];
                const auto phi = priVars[MassIndices::phaseFieldIdx];
                const auto cv = priVars[MassIndices::vaporIdx];
                const auto volume = scv.volume();

                liquidVolume += massProblem->liquidVolumeFraction(phi)*volume;
                evaporationMass += massProblem->evaporationMassSource(phi, cv)*volume;

                if (scv.dofPosition()[0] >= massProblem->interfacePositionAt(scv.dofPosition(), time))
                {
                    const auto cvExact = massProblem->analyticalVaporConcentration(scv.dofPosition(), time);
                    const auto err = cv - cvExact;
                    vaporErr2 += err*err*volume;
                    vaporNorm2 += cvExact*cvExact*volume;
                }
            }
        }

        Scalar crossSectionArea = 1.0;
        for (int dir = 1; dir < MassGridGeometry::GridView::dimensionworld; ++dir)
            crossSectionArea *= massGridGeometry->bBoxMax()[dir] - massGridGeometry->bBoxMin()[dir];

        const auto interfacePosition = massGridGeometry->bBoxMin()[0] + liquidVolume/crossSectionArea;
        qoi.gasLength = massProblem->openingPosition() - interfacePosition;
        qoi.analyticalGasLength = massProblem->analyticalGasLength(time);
        qoi.gasLengthError = qoi.gasLength - qoi.analyticalGasLength;
        qoi.liquidMass = massProblem->liquidDensity()*liquidVolume;
        qoi.evaporationRate = evaporationMass/crossSectionArea;
        qoi.analyticalEvaporationRate = massProblem->analyticalMassFlux(time);
        qoi.vaporL2Error = vaporNorm2 > 0.0 ? std::sqrt(vaporErr2/vaporNorm2) : 0.0;
        return qoi;
    };

    std::ofstream qoiFile(massProblem->name() + "_qoi.csv");
    qoiFile << "time,ell,ell_tube,ell_error,ell_time_l2,evap_rate,evap_rate_tube,liquid_mass,vapor_l2_rel\n";

    Scalar ellErr2TimeIntegral = 0.0;
    Scalar timeIntegral = 0.0;
    auto writeQoI = [&](const Scalar time, const Scalar elapsedDt)
    {
        const auto qoi = computeQoI(time);
        if (elapsedDt > 0.0)
        {
            ellErr2TimeIntegral += qoi.gasLengthError*qoi.gasLengthError*elapsedDt;
            timeIntegral += elapsedDt;
        }
        const auto ellL2 = timeIntegral > 0.0 ? std::sqrt(ellErr2TimeIntegral/timeIntegral) : std::abs(qoi.gasLengthError);

        qoiFile << time << ','
                << qoi.gasLength << ','
                << qoi.analyticalGasLength << ','
                << qoi.gasLengthError << ','
                << ellL2 << ','
                << qoi.evaporationRate << ','
                << qoi.analyticalEvaporationRate << ','
                << qoi.liquidMass << ','
                << qoi.vaporL2Error << '\n';

        std::cout << "\033[1;36m[stefan-tube] t=" << time
                  << " ell=" << qoi.gasLength
                  << " ell_tube=" << qoi.analyticalGasLength
                  << " err=" << qoi.gasLengthError
                  << " ell_L2(t)=" << ellL2
                  << " evap=" << qoi.evaporationRate
                  << " evap_tube=" << qoi.analyticalEvaporationRate
                  << " cv_L2rel=" << qoi.vaporL2Error
                  << "\033[0m" << std::endl;
    };
    writeQoI(0.0, 0.0);

    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager, timeLoop, xOld);

    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    const auto outputInterval = getParam<Scalar>("TimeLoop.OutputInterval", tEnd);
    Scalar nextOutputTime = outputInterval;

    timeLoop->start(); do
    {
        const auto elapsedDt = timeLoop->timeStepSize();
        nonLinearSolver.solve(x, *timeLoop);

        xOld = x;
        momentumGridVariables->advanceTimeStep();
        massGridVariables->advanceTimeStep();

        timeLoop->advanceTimeStep();
        writeQoI(timeLoop->time(), elapsedDt);

        if (timeLoop->time() >= nextOutputTime - 1e-10 || timeLoop->finished())
        {
            vtkWriter.write(timeLoop->time());
            while (timeLoop->time() >= nextOutputTime - 1e-10)
                nextOutputTime += outputInterval;
        }

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        static const int adaptEvery = getParam<int>("Adaptive.AdaptEveryNSteps", 1);
        if (getParam<bool>("Adaptive.EnableInTimeStep", true)
            && (timeLoop->timeStepIndex() % adaptEvery == 0))
        {
            indicator.calculate(x[massIdx], refineTol, coarsenTol);
            if (getParam<bool>("Adaptive.RefineAtVelocityGradient", true))
                indicator.includeVelocityJumps(
                    x[momentumIdx],
                    *momentumGridGeometry,
                    getParam<Scalar>("Adaptive.VelocityGradientRefineTolerance", refineTol),
                    getParam<Scalar>("Adaptive.VelocityGradientCoarsenTolerance", coarsenTol)
                );
            if (markElements(gridManager.grid(), indicator)
                && adapt(gridManager.grid(), dataTransfer))
            {
                xOld = x;
                couplingManager->updateSolution(x);
                massGridVariables->updateAfterGridAdaption(x[massIdx]);
                momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);
                couplingManager->init(momentumProblem, massProblem,
                                      std::make_tuple(momentumGridVariables, massGridVariables), x);
                assembler->updateAfterGridAdaption();
                std::cout << "\033[1;34m[adapt] step " << timeLoop->timeStepIndex()
                          << ": leafCells=" << gridManager.grid().leafGridView().size(0)
                          << "\033[0m" << std::endl;
            }
        }
    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << leafGridView.comm().size() << " processes." << std::endl;

    return 0;
}
