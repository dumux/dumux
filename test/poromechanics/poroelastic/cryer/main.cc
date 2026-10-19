// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Cryer sphere consolidation benchmark for linear poroelasticity.
 *
 * A porous sphere is loaded by a pressure on its drained surface. The pore pressure at the centre
 * first rises above its value right after loading and then decays (Mandel–Cryer effect).
 * The centre pressure is compared with Cryer's solution.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/experimental/assembly/multistagemultidomainfvassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/geometry/intersectingentities.hh>

#include <test/poromechanics/cryer_analytical.hh>

#include "properties.hh"
#include "volumetricstrain.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using OnePTypeTag = Properties::TTag::CryerOneP;
    using PoroMechTypeTag = Properties::TTag::CryerPoroElastic;
    using Scalar = GetPropType<OnePTypeTag, Properties::Scalar>;

    using GridManager = Dumux::GridManager<GetPropType<OnePTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using OnePGridDiscretization = GetPropType<OnePTypeTag, Properties::GridGeometry>;
    using PoroMechGridDiscretization = GetPropType<PoroMechTypeTag, Properties::GridGeometry>;
    auto onePGridDiscretization = std::make_shared<OnePGridDiscretization>(leafGridView);
    auto poroMechGridDiscretization = std::make_shared<PoroMechGridDiscretization>(leafGridView);

    using Traits = MultiDomainTraits<OnePTypeTag, PoroMechTypeTag>;
    using CouplingManager = GetPropType<OnePTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    using OnePProblem = GetPropType<OnePTypeTag, Properties::Problem>;
    using PoroMechProblem = GetPropType<PoroMechTypeTag, Properties::Problem>;
    auto onePSpatialParams = std::make_shared<typename OnePProblem::SpatialParams>(
        onePGridDiscretization, couplingManager);
    auto onePProblem = std::make_shared<OnePProblem>(onePGridDiscretization, onePSpatialParams, "OneP");
    auto poroMechSpatialParams = std::make_shared<typename PoroMechProblem::SpatialParams>(
        poroMechGridDiscretization, couplingManager);
    auto poroMechProblem = std::make_shared<PoroMechProblem>(
        poroMechGridDiscretization, poroMechSpatialParams, "PoroElastic");

    using SolutionVector = typename Traits::SolutionVector;
    static const auto onePId = Traits::template SubDomain<0>::Index();
    static const auto poroMechId = Traits::template SubDomain<1>::Index();
    SolutionVector x;
    x[onePId].resize(onePGridDiscretization->numDofs());
    x[poroMechId].resize(poroMechGridDiscretization->numDofs());
    onePProblem->applyInitialSolution(x[onePId]);
    poroMechProblem->applyInitialSolution(x[poroMechId]);
    SolutionVector xOld = x;

    couplingManager->init(onePProblem, poroMechProblem, x);

    using OnePGridVariables = GetPropType<OnePTypeTag, Properties::GridVariables>;
    using PoroMechGridVariables = GetPropType<PoroMechTypeTag, Properties::GridVariables>;
    auto onePGridVariables = std::make_shared<OnePGridVariables>(onePProblem, onePGridDiscretization);
    auto poroMechGridVariables = std::make_shared<PoroMechGridVariables>(poroMechProblem, poroMechGridDiscretization);
    onePGridVariables->init(x[onePId]);
    poroMechGridVariables->init(x[poroMechId]);

    const Scalar K = getParam<Scalar>("SpatialParams.BulkModulus");
    const Scalar G = getParam<Scalar>("SpatialParams.ShearModulus");
    const Scalar alphaB = getParam<Scalar>("SpatialParams.BiotCoefficient", 1.0);
    const Scalar Sp = getParam<Scalar>("SpatialParams.StorageCoefficient", 0.0);
    const Scalar kappa = getParam<Scalar>("SpatialParams.InitialPermeability");
    const Scalar muf = getParam<Scalar>("Component.LiquidDynamicViscosity");
    const Scalar R0 = getParam<Scalar>("Grid.SphereRadius", 1.0);
    const Scalar p0 = getParam<Scalar>("Problem.PressureLoad");
    CryerAnalyticalSolution<Scalar> analytical(K, G, alphaB, Sp, kappa, muf, R0);

    std::vector<Scalar> pAnalytical(leafGridView.size(0), 0.0);
    std::vector<Scalar> divU(leafGridView.size(0), 0.0);
    const auto updateOutputFields = [&](Scalar t)
    {
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = onePGridDiscretization->elementMapper().index(element);
            const auto center = element.geometry().center();
            pAnalytical[eIdx] = t > 0.0 ? analytical.normalizedPressure(center.two_norm(), t)*p0 : 0.0;
            const auto elemSol = elementSolution(element, x[poroMechId], *poroMechGridDiscretization);
            divU[eIdx] = volumetricStrain(element, *poroMechGridDiscretization, elemSol, center);
        }
    };

    using OnePSolutionVector = GetPropType<OnePTypeTag, Properties::SolutionVector>;
    using PoroMechSolutionVector = GetPropType<PoroMechTypeTag, Properties::SolutionVector>;
    using OnePVtkOutputModule = VtkOutputModule<OnePGridVariables, OnePSolutionVector>;
    using PoroMechVtkOutputModule = VtkOutputModule<PoroMechGridVariables, PoroMechSolutionVector>;
    OnePVtkOutputModule onePVtkWriter(*onePGridVariables, x[onePId], onePProblem->name());
    PoroMechVtkOutputModule poroMechVtkWriter(*poroMechGridVariables, x[poroMechId], poroMechProblem->name());
    GetPropType<OnePTypeTag, Properties::IOFields>::initOutputModule(onePVtkWriter);
    GetPropType<PoroMechTypeTag, Properties::IOFields>::initOutputModule(poroMechVtkWriter);
    onePVtkWriter.addField(pAnalytical, "p_analytical", Vtk::FieldType::element);
    onePVtkWriter.addField(divU, "divU", Vtk::FieldType::element);

    updateOutputFields(0.0);
    onePVtkWriter.write(0.0);
    poroMechVtkWriter.write(0.0);

    const auto timeSteppingMethod = []() -> std::shared_ptr<Experimental::MultiStageMethod<Scalar>>
    {
        const auto scheme = getParam<std::string>("TimeLoop.Scheme", "ImplicitEuler");
        if (scheme == "ImplicitEuler")
            return std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
        else if (scheme == "CrankNicolson")
            return std::make_shared<Experimental::MultiStage::Theta<Scalar>>(0.5);
        else if (scheme == "DIRK3")
            return std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>();
        else
            DUNE_THROW(ParameterException, "Unknown TimeLoop.Scheme: " << scheme);
    }();

    using Assembler = Experimental::MultiStageMultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(onePProblem, poroMechProblem),
        std::make_tuple(onePGridDiscretization, poroMechGridDiscretization),
        std::make_tuple(onePGridVariables, poroMechGridVariables),
        couplingManager, timeSteppingMethod, xOld
    );

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    Experimental::MultiStageTimeStepper<NewtonSolver> timeStepper(nonLinearSolver, timeSteppingMethod);

    // the cell that contains the centre and whose centre is closest to it
    const typename OnePGridDiscretization::GlobalCoordinate origin(0.0);
    std::size_t centerElementIdx = 0;
    Scalar centerCellRadius = std::numeric_limits<Scalar>::max();
    for (const auto eIdx : intersectingEntities(origin, onePGridDiscretization->boundingBoxTree()))
    {
        const Scalar r = onePGridDiscretization->element(eIdx).geometry().center().two_norm();
        if (r < centerCellRadius)
        {
            centerCellRadius = r;
            centerElementIdx = eIdx;
        }
    }
    const auto centerElement = onePGridDiscretization->element(centerElementIdx);

    const Scalar tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.Dt");
    const Scalar dtMax = getParam<Scalar>("TimeLoop.MaxTimeStepSize", tEnd/100.0);
    const Scalar dtGrowth = getParam<Scalar>("TimeLoop.DtGrowthFactor", 1.0);

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(dtMax);

    std::ofstream csvOut("cryer_center_pressure.csv");
    csvOut << "# time[s], t*cv/R0^2, p_f_center[Pa], p_analytical[Pa], p_f/p0, p_analytical/p0, J_center\n";

    Scalar maxCenterDeviation = 0.0;

    timeLoop->start();
    do
    {
        timeStepper.step(x, timeLoop->time(), timeLoop->timeStepSize());
        xOld = x;

        timeLoop->advanceTimeStep();
        onePGridVariables->advanceTimeStep();
        poroMechGridVariables->advanceTimeStep();

        const Scalar t = timeLoop->time();
        const Scalar pCenter = x[onePId][centerElementIdx][0];
        const Scalar pAnalyticalCenter = analytical.normalizedPressure(centerCellRadius, t)*p0;
        const auto centerElemSol = elementSolution(centerElement, x[poroMechId], *poroMechGridDiscretization);
        const Scalar JCenter = 1.0 + volumetricStrain(centerElement, *poroMechGridDiscretization,
                                                      centerElemSol, centerElement.geometry().center());
        maxCenterDeviation = std::max(maxCenterDeviation, std::abs(pCenter - pAnalyticalCenter)/p0);

        csvOut << t << ", " << t*analytical.cv()/(R0*R0)
               << ", " << pCenter << ", " << pAnalyticalCenter
               << ", " << pCenter/p0 << ", " << pAnalyticalCenter/p0
               << ", " << JCenter << "\n";

        updateOutputFields(t);
        onePVtkWriter.write(t);
        poroMechVtkWriter.write(t);
        timeLoop->reportTimeStep();

        dt = std::min(dt*dtGrowth, dtMax);
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
