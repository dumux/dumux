// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Case A planar Stefan evaporation verification draft.
 */

#include <config.h>

#include <cmath>
#include <fstream>
#include <iostream>

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
struct CaseAQoI
{
    Scalar interfacePosition = 0.0;
    Scalar analyticalInterfacePosition = 0.0;
    Scalar interfacePositionError = 0.0;
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

    const auto initMaxLevel = getParam<std::size_t>("Adaptive.InitMaxLevel", 0);
    for (std::size_t i = 0; i < initMaxLevel; ++i)
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
    }

    // The initial AMR loop transfers a profile that was sampled on the coarser
    // previous grid. Re-sample the analytical initial state on the final initial
    // grid so the Stefan convergence study is not polluted by transfer error.
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    xOld = x;

    couplingManager->init(momentumProblem, massProblem,
                          std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->updateAfterGridAdaption(x[massIdx]);
    momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);

    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    auto computeQoI = [&](const Scalar time)
    {
        CaseAQoI<Scalar> qoi;

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

                if (scv.dofPosition()[0] >= massProblem->analyticalInterfacePosition(time))
                {
                    const auto cvExact = massProblem->analyticalVaporConcentration(scv.dofPosition(), time);
                    const auto err = cv - cvExact;
                    vaporErr2 += err*err*volume;
                    vaporNorm2 += cvExact*cvExact*volume;
                }
            }
        }

        const auto height = massGridGeometry->bBoxMax()[1] - massGridGeometry->bBoxMin()[1];
        qoi.interfacePosition = massGridGeometry->bBoxMin()[0] + liquidVolume/height;
        qoi.analyticalInterfacePosition = massProblem->analyticalInterfacePosition(time);
        qoi.interfacePositionError = qoi.interfacePosition - qoi.analyticalInterfacePosition;
        qoi.liquidMass = massProblem->liquidDensity()*liquidVolume;
        qoi.evaporationRate = evaporationMass/height;
        qoi.analyticalEvaporationRate = massProblem->analyticalMassFlux(time);
        qoi.vaporL2Error = vaporNorm2 > 0.0 ? std::sqrt(vaporErr2/vaporNorm2) : 0.0;
        return qoi;
    };

    std::ofstream qoiFile(massProblem->name() + "_qoi.csv");
    qoiFile << "time,xi,xi_stefan,xi_error,xi_time_l2,evap_rate,evap_rate_stefan,liquid_mass,vapor_l2_rel\n";

    Scalar xiErr2TimeIntegral = 0.0;
    Scalar timeIntegral = 0.0;
    auto writeQoI = [&](const Scalar time, const Scalar elapsedDt)
    {
        const auto qoi = computeQoI(time);
        if (elapsedDt > 0.0)
        {
            xiErr2TimeIntegral += qoi.interfacePositionError*qoi.interfacePositionError*elapsedDt;
            timeIntegral += elapsedDt;
        }
        const auto xiL2 = timeIntegral > 0.0 ? std::sqrt(xiErr2TimeIntegral/timeIntegral) : std::abs(qoi.interfacePositionError);

        qoiFile << time << ','
                << qoi.interfacePosition << ','
                << qoi.analyticalInterfacePosition << ','
                << qoi.interfacePositionError << ','
                << xiL2 << ','
                << qoi.evaporationRate << ','
                << qoi.analyticalEvaporationRate << ','
                << qoi.liquidMass << ','
                << qoi.vaporL2Error << '\n';

        std::cout << "\033[1;36m[stefan] t=" << time
                  << " xi=" << qoi.interfacePosition
                  << " xi_stefan=" << qoi.analyticalInterfacePosition
                  << " err=" << qoi.interfacePositionError
                  << " xi_L2(t)=" << xiL2
                  << " evap=" << qoi.evaporationRate
                  << " evap_stefan=" << qoi.analyticalEvaporationRate
                  << " cv_L2rel=" << qoi.vaporL2Error
                  << "\033[0m" << std::endl;
    };
    writeQoI(0.0, 0.0);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager, timeLoop, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
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
