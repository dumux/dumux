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
 * A porous sphere of radius R₀ = 1 m is subjected to a uniform surface pressure p₀.
 * Due to the Mandel-Cryer effect, the pore pressure at the centre initially rises
 * ABOVE p₀ before draining to zero.
 *
 * The benchmark runs three pressure loads (p₀/K = 0.0001, 0.25, 0.5) and two
 * permeability models (constant / Kozeny–Carman) and compares against the small-
 * deformation analytical solution of Cryer (1963).
 *
 * Usage (for one load, via command line or params.input):
 *   ./test_cryer -Problem.PressureLoad 100.0 -SpatialParams.PermeabilityModel KozenyCarman
 */
#include <config.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
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

#include "properties.hh"
#include "cryer_analytical.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomTT  = Properties::TTag::CryerMomentum;
    using SpTT   = Properties::TTag::CryerSolidPressure;
    using FpTT   = Properties::TTag::CryerFluidPressure;
    using MDTraits = Properties::CryerMDTraits;
    using Scalar = typename MDTraits::Scalar;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    // -----------------------------------------------------------------------
    // Grid
    // -----------------------------------------------------------------------
    using Grid = GetPropType<MomTT, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    // -----------------------------------------------------------------------
    // Grid geometries
    // -----------------------------------------------------------------------
    using MomGG = GetPropType<MomTT, Properties::GridGeometry>;
    using SpGG  = GetPropType<SpTT,  Properties::GridGeometry>;
    using FpGG  = GetPropType<FpTT,  Properties::GridGeometry>;
    auto momGG = std::make_shared<MomGG>(leafGridView);
    auto spGG  = std::make_shared<SpGG>(leafGridView);
    auto fpGG  = std::make_shared<FpGG>(leafGridView);

    // -----------------------------------------------------------------------
    // Coupling manager
    // -----------------------------------------------------------------------
    using CouplingManager = PoroElasticLargeDefCouplingManager<MDTraits>;
    auto couplingManager = std::make_shared<CouplingManager>(momGG, spGG, fpGG);

    // -----------------------------------------------------------------------
    // Problems (spatial params are created internally by FVProblemWithSpatialParams)
    // -----------------------------------------------------------------------
    using MomProblem = GetPropType<MomTT, Properties::Problem>;
    using SpProblem  = GetPropType<SpTT,  Properties::Problem>;
    using FpProblem  = GetPropType<FpTT,  Properties::Problem>;
    auto momProblem = std::make_shared<MomProblem>(momGG, couplingManager);
    auto spProblem  = std::make_shared<SpProblem>(spGG,   couplingManager);
    auto fpProblem  = std::make_shared<FpProblem>(fpGG,   couplingManager);

    // -----------------------------------------------------------------------
    // Solution vector
    // -----------------------------------------------------------------------
    using SolVec = typename MDTraits::SolutionVector;
    SolVec x;
    x[CouplingManager::momentumIdx].resize(momGG->numDofs());
    x[CouplingManager::solidPressureIdx].resize(spGG->numDofs());
    x[CouplingManager::fluidPressureIdx].resize(fpGG->numDofs());
    x = 0.0;
    SolVec xOld = x;

    // -----------------------------------------------------------------------
    // Grid variables
    // -----------------------------------------------------------------------
    using MomGV = GetPropType<MomTT, Properties::GridVariables>;
    using SpGV  = GetPropType<SpTT,  Properties::GridVariables>;
    using FpGV  = GetPropType<FpTT,  Properties::GridVariables>;
    auto momGV = std::make_shared<MomGV>(momProblem, momGG);
    auto spGV  = std::make_shared<SpGV>(spProblem,   spGG);
    auto fpGV  = std::make_shared<FpGV>(fpProblem,   fpGG);
    momGV->init(x[CouplingManager::momentumIdx]);
    spGV->init(x[CouplingManager::solidPressureIdx]);
    fpGV->init(x[CouplingManager::fluidPressureIdx]);

    couplingManager->init(momProblem, spProblem, fpProblem, x);
    couplingManager->computeColorsForAssembly();

    // -----------------------------------------------------------------------
    // Time stepping method
    // -----------------------------------------------------------------------
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

    // -----------------------------------------------------------------------
    // Assembler + linear/Newton solvers
    // -----------------------------------------------------------------------
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

    // -----------------------------------------------------------------------
    // Analytical solution
    // -----------------------------------------------------------------------
    const Scalar K = momProblem->spatialParams().bulkModulus();
    const Scalar G = momProblem->spatialParams().shearModulus();
    const Scalar kappa0 = getParam<Scalar>("SpatialParams.InitialPermeability");
    const Scalar muf    = getParam<Scalar>("SpatialParams.FluidViscosity");
    const Scalar alphaB = getParam<Scalar>("SpatialParams.BiotCoefficient", 1.0);
    const Scalar Sp     = getParam<Scalar>("SpatialParams.StorageCoefficient", 0.0);
    const Scalar R0     = getParam<Scalar>("Grid.SphereRadius", 1.0);
    const Scalar p0     = getParam<Scalar>("Problem.PressureLoad");

    CryerAnalyticalSolution<Scalar> analytical(K, G, alphaB, Sp, kappa0, muf, R0);

    // -----------------------------------------------------------------------
    // VTK output
    // -----------------------------------------------------------------------
    VtkOutputModule<FpGV, GetPropType<FpTT, Properties::SolutionVector>>
        vtkWriter(*fpGV, x[CouplingManager::fluidPressureIdx], "cryer_fluid_pressure");
    vtkWriter.addVolumeVariable([](const auto& v){ return v.fluidPressure(); }, "p_f");

    // Analytical pressure field p(r,t)/p0 as VTK field
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

    // -----------------------------------------------------------------------
    // Time loop
    // -----------------------------------------------------------------------
    const Scalar tEnd   = getParam<Scalar>("TimeLoop.TEnd");
    auto dt             = getParam<Scalar>("TimeLoop.Dt");
    const Scalar dtMax  = getParam<Scalar>("TimeLoop.MaxTimeStepSize", tEnd / 100.0);
    // Geometric time-step growth (>=1). With dtGrowth > 1 the time steps are
    // log-spaced (small dt early to resolve the pressure rise, large dt later
    // for the slow decay), as anticipated in params.input. Default 1.0 keeps the
    // fixed-dt behaviour used by the CI smoke test.
    const Scalar dtGrowth = getParam<Scalar>("TimeLoop.DtGrowthFactor", 1.0);

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(dtMax);

    updateAnalytical(0.0);
    vtkWriter.write(0.0);

    // Time series CSV for the centre pressure
    std::ofstream csvOut("cryer_center_pressure.csv");
    csvOut << "# time[s], t*cv/R0^2, p_f_center[Pa], p_analytical[Pa], p_f/p0, p_analytical/p0, J_center\n";

    // Find sphere centre element
    const typename MomGG::GlobalCoordinate centerPos(0.0);

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

        // Evaluate numerical centre pressure and volume ratio J = V/V_0
        Scalar pfCenter = 0.0;
        Scalar Jcenter = 1.0;
        const auto& fpSol = x[CouplingManager::fluidPressureIdx];
        const auto centerElements = intersectingEntities(centerPos, fpGG->boundingBoxTree());
        if (!centerElements.empty())
        {
            const auto elem = fpGG->element(centerElements[0]);
            const auto elemSol = elementSolution(elem, fpSol, *fpGG);
            pfCenter = evalSolution(elem, elem.geometry(), *fpGG, elemSol, centerPos)[0];

            // J = det(F) at the centre from the (just-converged) displacement field
            couplingManager->updateSolution(x);
            auto fpFvGeom = localView(*fpGG);
            fpFvGeom.bindElement(elem);
            Jcenter = couplingManager->deformationGradientAtPoint(fpFvGeom, centerPos).determinant();
        }

        const Scalar pAnal = analytical.normalizedCenterPressure(t) * p0;

        csvOut << t << ", " << tNorm
               << ", " << pfCenter << ", " << pAnal
               << ", " << pfCenter/p0 << ", " << pAnal/p0
               << ", " << Jcenter << "\n";

        updateAnalytical(t);
        vtkWriter.write(t);
        timeLoop->reportTimeStep();

        // grow the time step (geometric, capped at dtMax) for log-spaced sampling
        dt = std::min(dt * dtGrowth, dtMax);
        timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());
    nonLinearSolver->report();

    if (leafGridView.comm().rank() == 0)
        Parameters::print();

    return 0;
}
