// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Cryer sphere consolidation benchmark for the large-deformation poroelastic model.
 *
 * A porous sphere is loaded by a pressure on its drained surface. The pore pressure at the centre
 * first rises above its value right after loading and then decays (Mandel–Cryer effect).
 * The centre pressure is compared with Cryer's small-deformation solution.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/multistagemultidomainassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include <test/poromechanics/cryer_analytical.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomTT = Properties::TTag::CryerMomentum;
    using SpTT = Properties::TTag::CryerSolidPressure;
    using FpTT = Properties::TTag::CryerFluidPressure;
    using MDTraits = Properties::CryerMDTraits;
    using Scalar = typename MDTraits::Scalar;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<MomTT, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using MomGG = GetPropType<MomTT, Properties::GridGeometry>;
    using SpGG = GetPropType<SpTT, Properties::GridGeometry>;
    using FpGG = GetPropType<FpTT, Properties::GridGeometry>;
    auto momGG = std::make_shared<MomGG>(leafGridView);
    auto spGG = std::make_shared<SpGG>(leafGridView);
    auto fpGG = std::make_shared<FpGG>(leafGridView);

    using CouplingManager = PoroElasticLargeDefCouplingManager<MDTraits>;
    auto couplingManager = std::make_shared<CouplingManager>(momGG, spGG, fpGG);

    using MomProblem = GetPropType<MomTT, Properties::Problem>;
    using SpProblem = GetPropType<SpTT, Properties::Problem>;
    using FpProblem = GetPropType<FpTT, Properties::Problem>;
    auto momProblem = std::make_shared<MomProblem>(momGG, couplingManager);
    auto spProblem = std::make_shared<SpProblem>(spGG, couplingManager);
    auto fpProblem = std::make_shared<FpProblem>(fpGG, couplingManager);

    using SolVec = typename MDTraits::SolutionVector;
    SolVec x;
    x[CouplingManager::momentumIdx].resize(momGG->numDofs());
    x[CouplingManager::solidPressureIdx].resize(spGG->numDofs());
    x[CouplingManager::fluidPressureIdx].resize(fpGG->numDofs());
    x = 0.0;
    SolVec xOld = x;

    using MomGV = GetPropType<MomTT, Properties::GridVariables>;
    using SpGV = GetPropType<SpTT, Properties::GridVariables>;
    using FpGV = GetPropType<FpTT, Properties::GridVariables>;
    auto momGV = std::make_shared<MomGV>(momProblem, momGG);
    auto spGV = std::make_shared<SpGV>(spProblem, spGG);
    auto fpGV = std::make_shared<FpGV>(fpProblem, fpGG);
    momGV->init(x[CouplingManager::momentumIdx]);
    spGV->init(x[CouplingManager::solidPressureIdx]);
    fpGV->init(x[CouplingManager::fluidPressureIdx]);

    couplingManager->init(momProblem, spProblem, fpProblem, x);
    couplingManager->computeColorsForAssembly();

    const auto tsScheme = getParam<std::string>("TimeLoop.Scheme", "ImplicitEuler");
    std::shared_ptr<Experimental::MultiStageMethod<Scalar>> msMethod;
    if (tsScheme == "ImplicitEuler")
        msMethod = std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
    else if (tsScheme == "CrankNicolson")
        msMethod = std::make_shared<Experimental::MultiStage::Theta<Scalar>>(0.5);
    else if (tsScheme == "DIRK3")
        msMethod = std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>();
    else
        DUNE_THROW(ParameterException, "Unknown TimeLoop.Scheme: " << tsScheme);

    using Assembler = Experimental::MultiStageMultiDomainAssembler<MDTraits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(momProblem, spProblem, fpProblem),
        std::make_tuple(momGG, spGG, fpGG),
        std::make_tuple(momGV, spGV, fpGV),
        couplingManager, msMethod, xOld
    );
    assembler->setLinearSystem();

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits,
                                           LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    using TimeStepper = Experimental::MultiStageTimeStepper<NewtonSolver>;
    TimeStepper timeStepper(nonLinearSolver, msMethod);

    const Scalar K = momProblem->spatialParams().bulkModulus();
    const Scalar G = momProblem->spatialParams().shearModulus();
    const Scalar kappa0 = getParam<Scalar>("SpatialParams.InitialPermeability");
    const Scalar muf = getParam<Scalar>("SpatialParams.FluidViscosity");
    const Scalar alphaB = getParam<Scalar>("SpatialParams.BiotCoefficient", 1.0);
    const Scalar Sp = getParam<Scalar>("SpatialParams.StorageCoefficient", 0.0);
    const Scalar R0 = getParam<Scalar>("Grid.SphereRadius", 1.0);
    const Scalar p0 = getParam<Scalar>("Problem.PressureLoad");

    CryerAnalyticalSolution<Scalar> analytical(K, G, alphaB, Sp, kappa0, muf, R0);

    VtkOutputModule<FpGV, GetPropType<FpTT, Properties::SolutionVector>>
        vtkWriter(*fpGV, x[CouplingManager::fluidPressureIdx], "cryer_fluid_pressure");
    vtkWriter.addVolumeVariable([](const auto& v){ return v.fluidPressure(); }, "p_f");

    std::vector<Scalar> pAnalytical(fpGG->gridView().size(0), 0.0);
    auto updateAnalytical = [&](Scalar t)
    {
        for (const auto& element : elements(fpGG->gridView()))
        {
            const auto eIdx = fpGG->elementMapper().index(element);
            const Scalar r = element.geometry().center().two_norm();
            pAnalytical[eIdx] = analytical.normalizedPressure(r, t) * p0;
        }
    };
    vtkWriter.addField(pAnalytical, "p_f_analytical", Dumux::Vtk::FieldType::element);

    const Scalar tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.Dt");
    const Scalar dtMax = getParam<Scalar>("TimeLoop.MaxTimeStepSize", tEnd / 100.0);
    const Scalar dtGrowth = getParam<Scalar>("TimeLoop.DtGrowthFactor", 1.0);

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(dtMax);

    updateAnalytical(0.0);
    vtkWriter.write(0.0);

    std::ofstream csvOut("cryer_center_pressure.csv");
    csvOut << "# time[s], t*cv/R0^2, p_f_center[Pa], p_analytical[Pa], p_f/p0, p_analytical/p0, J_center\n";

    const typename MomGG::GlobalCoordinate centerPos(0.0);
    Scalar maxCenterDeviation = 0.0;

    timeLoop->start();
    do
    {
        xOld = x;
        assembler->setPreviousSolution(xOld);

        timeStepper.step(x, timeLoop->time(), timeLoop->timeStepSize());

        timeLoop->advanceTimeStep();
        momGV->advanceTimeStep();
        spGV->advanceTimeStep();
        fpGV->advanceTimeStep();

        const Scalar t = timeLoop->time();
        const Scalar tNorm = t * analytical.cv() / (R0 * R0);

        Scalar pfCenter = 0.0;
        Scalar Jcenter = 1.0;
        const auto& fpSol = x[CouplingManager::fluidPressureIdx];
        const auto centerElements = intersectingEntities(centerPos, fpGG->boundingBoxTree());
        if (!centerElements.empty())
        {
            const auto elem = fpGG->element(centerElements[0]);
            const auto elemSol = elementSolution(elem, fpSol, *fpGG);
            pfCenter = evalSolution(elem, elem.geometry(), *fpGG, elemSol, centerPos)[0];

            couplingManager->updateSolution(x);
            auto fpElemDisc = localView(*fpGG);
            fpElemDisc.bindElement(elem);
            const CVFE::InterpolationPointData centerIpData(elem.geometry().local(centerPos), centerPos);
            Jcenter = couplingManager->deformationGradient(fpElemDisc, centerIpData).determinant();
        }

        const Scalar pAnal = analytical.normalizedCenterPressure(t) * p0;
        maxCenterDeviation = std::max(maxCenterDeviation, std::abs(pfCenter - pAnal)/p0);

        csvOut << t << ", " << tNorm
               << ", " << pfCenter << ", " << pAnal
               << ", " << pfCenter/p0 << ", " << pAnal/p0
               << ", " << Jcenter << "\n";

        updateAnalytical(t);
        vtkWriter.write(t);
        timeLoop->reportTimeStep();

        dt = std::min(dt * dtGrowth, dtMax);
        timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());
    nonLinearSolver->report();

    std::cout << "Maximum deviation of the centre pressure from Cryer's solution: "
              << maxCenterDeviation << " p0" << std::endl;
    if (hasParam("Problem.CenterPressureTolerance"))
    {
        const auto tolerance = getParam<Scalar>("Problem.CenterPressureTolerance");
        if (!(maxCenterDeviation <= tolerance))
            DUNE_THROW(Dune::Exception, "Centre pressure deviates from Cryer's solution by "
                        << maxCenterDeviation << " p0 > " << tolerance << " p0");
    }

    if (leafGridView.comm().rank() == 0)
        Parameters::print();

    return 0;
}
