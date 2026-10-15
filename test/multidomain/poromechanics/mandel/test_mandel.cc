// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Mandel's problem: consolidation of a poroelastic slab between two rigid plates,
 *        solved with the coupled single-phase flow and poroelastic models.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/experimental/assembly/multistagemultidomainfvassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>
#include <dumux/experimental/timestepping/generalizedalpha.hh>

#include <test/experimental/timestepping/implicitsolvecheck.hh>

#include "properties.hh"
#include "volumetricstrain.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using OnePTypeTag = Properties::TTag::OnePSub;
    using PoroMechTypeTag = Properties::TTag::PoroElasticSub;
    using Scalar = GetPropType<OnePTypeTag, Properties::Scalar>;

    using GridManager = Dumux::GridManager<GetPropType<OnePTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using OnePGridGeometry = GetPropType<OnePTypeTag, Properties::GridGeometry>;
    using PoroMechGridGeometry = GetPropType<PoroMechTypeTag, Properties::GridGeometry>;
    auto onePGridGeometry = std::make_shared<OnePGridGeometry>(leafGridView);
    auto poroMechGridGeometry = std::make_shared<PoroMechGridGeometry>(leafGridView);

    using Traits = MultiDomainTraits<OnePTypeTag, PoroMechTypeTag>;
    using CouplingManager = PoroMechanicsCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    auto analyticalSolution = std::make_shared<const MandelAnalyticalSolution<Scalar>>();

    using OnePProblem = GetPropType<OnePTypeTag, Properties::Problem>;
    using PoroMechProblem = GetPropType<PoroMechTypeTag, Properties::Problem>;
    auto onePSpatialParams = std::make_shared<typename OnePProblem::SpatialParams>(onePGridGeometry, couplingManager, analyticalSolution);
    auto onePProblem = std::make_shared<OnePProblem>(onePGridGeometry, onePSpatialParams, "OneP");
    auto poroMechSpatialParams = std::make_shared<typename PoroMechProblem::SpatialParams>(poroMechGridGeometry, couplingManager, analyticalSolution);
    auto poroMechProblem = std::make_shared<PoroMechProblem>(poroMechGridGeometry, poroMechSpatialParams, "PoroElastic");

    using SolutionVector = typename Traits::SolutionVector;
    static const auto onePId = Traits::template SubDomain<0>::Index();
    static const auto poroMechId = Traits::template SubDomain<1>::Index();
    SolutionVector x;
    x[onePId].resize(onePGridGeometry->numDofs());
    x[poroMechId].resize(poroMechGridGeometry->numDofs());
    onePProblem->applyInitialSolution(x[onePId]);
    poroMechProblem->applyInitialSolution(x[poroMechId]);
    SolutionVector xOld = x;

    couplingManager->init(onePProblem, poroMechProblem, x);

    using OnePGridVariables = GetPropType<OnePTypeTag, Properties::GridVariables>;
    using PoroMechGridVariables = GetPropType<PoroMechTypeTag, Properties::GridVariables>;
    auto onePGridVariables = std::make_shared<OnePGridVariables>(onePProblem, onePGridGeometry);
    auto poroMechGridVariables = std::make_shared<PoroMechGridVariables>(poroMechProblem, poroMechGridGeometry);
    onePGridVariables->init(x[onePId]);
    poroMechGridVariables->init(x[poroMechId]);

    auto pExact = x[onePId];
    auto uExact = x[poroMechId];
    std::vector<Scalar> divU(leafGridView.size(0));
    const auto updateOutputFields = [&](const Scalar time)
    {
        Dumux::parallelFor(leafGridView.size(0), [&](const std::size_t eIdx)
        {
            const auto element = onePGridGeometry->element(eIdx);
            const auto center = element.geometry().center();
            pExact[eIdx] = time > 0.0 ? analyticalSolution->pressure(center, time)
                                      : analyticalSolution->initialPressure(center);
            const auto elemSol = elementSolution(element, x[poroMechId], *poroMechGridGeometry);
            divU[eIdx] = volumetricStrain(element, *poroMechGridGeometry, elemSol, center);
        });

        for (const auto& vertex : vertices(leafGridView))
        {
            const auto vIdx = poroMechGridGeometry->vertexMapper().index(vertex);
            const auto pos = vertex.geometry().center();
            uExact[vIdx] = time > 0.0 ? analyticalSolution->displacement(pos, time)
                                      : analyticalSolution->initialDisplacement(pos);
        }
    };

    using OnePVtkOutputModule = VtkOutputModule<OnePGridVariables, GetPropType<OnePTypeTag, Properties::SolutionVector>>;
    using PoroMechVtkOutputModule = VtkOutputModule<PoroMechGridVariables, GetPropType<PoroMechTypeTag, Properties::SolutionVector>>;
    OnePVtkOutputModule onePVtkWriter(*onePGridVariables, x[onePId], onePProblem->name());
    PoroMechVtkOutputModule poroMechVtkWriter(*poroMechGridVariables, x[poroMechId], poroMechProblem->name());
    GetPropType<OnePTypeTag, Properties::IOFields>::initOutputModule(onePVtkWriter);
    GetPropType<PoroMechTypeTag, Properties::IOFields>::initOutputModule(poroMechVtkWriter);
    onePVtkWriter.addVolumeVariable([](const auto& v){ return v.porosity(); }, "porosity");
    onePVtkWriter.addVolumeVariable([](const auto& v){ return v.density(); }, "rho");
    onePVtkWriter.addField(pExact, "pExact");
    onePVtkWriter.addField(divU, "divU", Vtk::FieldType::element);
    poroMechVtkWriter.addField(uExact, "uExact");

    updateOutputFields(0.0);
    onePVtkWriter.write(0.0);
    poroMechVtkWriter.write(0.0);

    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.Dt");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxTimeStepSize"));
    const auto outputTimes = getParam<std::vector<Scalar>>("TimeLoop.CheckPoints");
    timeLoop->setCheckPoint(outputTimes.begin(), outputTimes.end());

    using Assembler = Experimental::MultiStageMultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    constexpr bool withInertia = Assembler::hasInertia();

    // the elastic waves are damped with the Chung-Hulbert method, the consolidation is resolved
    const auto generalizedAlpha = []
    { return std::make_shared<Experimental::GeneralizedAlpha::ChungHulbert<Scalar>>(getParam<Scalar>("TimeLoop.SpectralRadius")); };

    const auto multiStageMethod = []() -> std::shared_ptr<Experimental::MultiStageMethod<Scalar>>
    {
        const auto scheme = getParam<std::string>("TimeLoop.Scheme", "ImplicitEuler");
        if (scheme == "ImplicitEuler")
            return std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
        else if (scheme == "CrankNicolson")
            return std::make_shared<Experimental::MultiStage::Theta<Scalar>>(0.5);
        else if (scheme == "DIRK3")
            return std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>();
        else
            DUNE_THROW(ParameterException, "Unknown time stepping scheme " << scheme
                       << ". Possible values are ImplicitEuler, CrankNicolson and DIRK3.");
    };

    const auto timeSteppingMethod = [&]
    {
        if constexpr (withInertia)
            return generalizedAlpha();
        else
            return multiStageMethod();
    }();

    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(onePProblem, poroMechProblem),
        std::make_tuple(onePGridGeometry, poroMechGridGeometry),
        std::make_tuple(onePGridVariables, poroMechGridVariables),
        couplingManager, timeSteppingMethod, xOld
    );

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    auto timeStepper = [&](const auto& method)
    {
        if constexpr (withInertia)
        {
            // the solid starts at rest
            Experimental::GeneralizedAlphaTimeStepper<NewtonSolver> stepper(nonLinearSolver, method);
            stepper.initialize(x, 0.0);
            return stepper;
        }
        else
            return Experimental::MultiStageTimeStepper<NewtonSolver>(nonLinearSolver, method);
    }(timeSteppingMethod);

    const Scalar inverseDt = 1.0/timeLoop->timeStepSize();
    if (Test::checkImplicitSolveContract(*assembler, x, timeLoop->time(), inverseDt, 1.0, withInertia ? inverseDt*inverseDt : 0.0) == 0)
        DUNE_THROW(Dune::Exception, "Expected Dirichlet constraints in the implicit solve contract check");

    // a translation of the solid changes neither strain nor porosity, so the inertia summed over
    // all degrees of freedom changes by the density of the mixture times the volume times the translation
    if constexpr (withInertia)
    {
        using Displacement = Dune::FieldVector<Scalar, 2>;
        const auto summedInertia = [&](const SolutionVector& sol)
        {
            couplingManager->updateSolution(sol);
            assembler->updateGridVariables(sol);
            typename Assembler::ResidualType inertia;
            assembler->assembleOperators(sol, Experimental::TimeLevel<Scalar>(0.0), {.inertia = &inertia});
            Displacement total(0.0);
            for (const auto& block : inertia[poroMechId])
                total += block;
            return total;
        };

        const Displacement translation{0.3, -0.2};
        auto xTranslated = x;
        for (auto& block : xTranslated[poroMechId])
            block += translation;

        auto change = summedInertia(xTranslated);
        change -= summedInertia(x);

        const auto phi = getParam<Scalar>("Problem.InitialPorosity");
        const auto mixtureDensity = phi*getParam<Scalar>("Component.LiquidDensity")
                                    + (1.0 - phi)*getParam<Scalar>("Component.SolidDensity");
        auto volume = poroMechGridGeometry->bBoxMax() - poroMechGridGeometry->bBoxMin();
        auto expected = translation;
        expected *= mixtureDensity*volume[0]*volume[1];
        auto difference = change;
        difference -= expected;
        std::cout << "[Inertia] translation: change of the summed inertia " << change << ", expected " << expected << std::endl;
        if (!(difference.infinity_norm() < 1e-10*expected.infinity_norm()))
            DUNE_THROW(Dune::Exception, "The inertia of a translation is not the mixture density times the volume");
    }

    timeLoop->start(); do
    {
        timeStepper.step(x, *timeLoop);
        xOld = x;

        timeLoop->advanceTimeStep();
        onePGridVariables->advanceTimeStep();
        poroMechGridVariables->advanceTimeStep();

        if (timeLoop->isCheckPoint() || timeLoop->finished())
        {
            updateOutputFields(timeLoop->time());
            onePVtkWriter.write(timeLoop->time());
            poroMechVtkWriter.write(timeLoop->time());
        }

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    {
        const Scalar t = timeLoop->time();
        Scalar errorP = 0.0;
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = onePGridGeometry->elementMapper().index(element);
            const auto delta = x[onePId][eIdx][0] - analyticalSolution->pressure(element.geometry().center(), t);
            errorP += element.geometry().volume()*delta*delta;
        }

        Scalar errorUx = 0.0, errorUy = 0.0;
        auto fvGeometry = localView(*poroMechGridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto u = analyticalSolution->displacement(scv.dofPosition(), t);
                const auto deltaX = x[poroMechId][scv.dofIndex()][0] - u[0];
                const auto deltaY = x[poroMechId][scv.dofIndex()][1] - u[1];
                errorUx += scv.volume()*deltaX*deltaX;
                errorUy += scv.volume()*deltaY*deltaY;
            }
        }

        using std::sqrt;
        errorP = sqrt(errorP);
        errorUx = sqrt(errorUx);
        errorUy = sqrt(errorUy);

        if (leafGridView.comm().rank() == 0)
        {
            std::cout << "** L2 errors at t = " << t
                      << ": L2(p) = " << errorP
                      << ", L2(ux) = " << errorUx
                      << ", L2(uy) = " << errorUy << std::endl;
        }
    }

    nonLinearSolver->report();
    timeLoop->finalize(leafGridView.comm());

    if (leafGridView.comm().rank() == 0)
        Parameters::print();

    return 0;
}
