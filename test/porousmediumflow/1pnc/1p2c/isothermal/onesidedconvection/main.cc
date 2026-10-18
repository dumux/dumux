// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief One-sided convection: dissolution flux of the diffusive boundary layer and amplification of a finger mode
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::OneSidedConvection;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    constexpr int soluteIdx = FluidSystem::comp1Idx;

    // Rayleigh number Ra = U H / D_eff with the buoyancy velocity U = K drho g / mu, time scale phi H / U;
    // the parameters are constant, so the volume variables of any solution provide them
    const Scalar height = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];
    const Scalar length = gridGeometry->bBoxMax()[0] - gridGeometry->bBoxMin()[0];
    auto fvGeometry = localView(*gridGeometry);
    auto elemVolVars = localView(gridVariables->curGridVolVars());
    const auto element = *elements(leafGridView).begin();
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, x);
    const auto& volVars = elemVolVars[*scvs(fvGeometry).begin()];
    const Scalar effectiveDiffusivity = volVars.effectiveDiffusionCoefficient(0, FluidSystem::comp0Idx, soluteIdx);
    const Scalar g = -problem->spatialParams().gravity(gridGeometry->bBoxMin())[1];
    const Scalar buoyancyVelocity = volVars.permeability()*problem->densityDifference()*g/volVars.viscosity(0);
    const Scalar rayleigh = buoyancyVelocity*height/effectiveDiffusivity;
    const Scalar timeScale = volVars.porosity()*height/buoyancyVelocity;

    problem->setRayleighNumber(rayleigh);
    problem->applyInitialSolution(x);
    auto xOld = x;
    gridVariables->init(x);

    const Scalar k = M_PI*height/length;
    const bool seeded = getParam<Scalar>("Problem.Amplitude") != 0.0;

    // amplitude of cos(kx) on each row of vertices of the tensor-product grid
    struct Row
    {
        Scalar volume = 0.0, mass = 0.0, cosine = 0.0, massCosine = 0.0, cosineSquared = 0.0;
        Scalar min = std::numeric_limits<Scalar>::max(), max = std::numeric_limits<Scalar>::lowest();
    };
    const auto rowsOf = [&](const SolutionVector& sol)
    {
        std::map<Scalar, Row> rows;
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& pos = scv.dofPosition();
                const Scalar c = std::cos(M_PI*(pos[0] - gridGeometry->bBoxMin()[0])/length);
                const Scalar massFraction = sol[scv.dofIndex()][soluteIdx];
                auto& row = rows[pos[1]];
                row.volume += scv.volume();
                row.mass += massFraction*scv.volume();
                row.cosine += c*scv.volume();
                row.massCosine += massFraction*c*scv.volume();
                row.cosineSquared += c*c*scv.volume();
                row.min = std::min(row.min, massFraction);
                row.max = std::max(row.max, massFraction);
            }
        }
        return rows;
    };

    const auto rowAmplitude = [](const Row& row)
    { return (row.massCosine - row.mass/row.volume*row.cosine)/row.cosineSquared; };

    // L2 norm of the mode amplitude over the depth, mean mass fraction, largest horizontal variation in a row
    const auto measure = [&](const SolutionVector& sol)
    {
        Scalar normSquared = 0.0, mass = 0.0, rowDeviation = 0.0;
        for (const auto& [z, row] : rowsOf(sol))
        {
            normSquared += rowAmplitude(row)*rowAmplitude(row)*row.volume/(length*height);
            mass += row.mass/(length*height);
            rowDeviation = std::max(rowDeviation, row.max - row.min);
        }
        return std::array<Scalar, 3>{std::sqrt(normSquared), mass, rowDeviation};
    };

    // mean mass fraction of the diffusive profile in the semi-infinite layer at dimensionless time t
    const auto exactMass = [&](Scalar t) { return 2.0*std::sqrt(t/(M_PI*rayleigh)); };

    const Scalar t0 = problem->initialTime();
    const Scalar dt = getParam<Scalar>("TimeLoop.DimensionlessTimeStepSize")*timeScale;
    const Scalar tEnd = (getParam<Scalar>("TimeLoop.DimensionlessEndTime") - t0)*timeScale;
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(dt);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    VtkOutputModule vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    std::ofstream logFile(problem->name() + ".csv");
    logFile << std::setprecision(12) << "time,amplitude,flux,exactFlux,maxRowDeviation\n";
    auto measured = measure(x);
    Scalar amplitude = measured[0], mass = measured[1], rowDeviation = measured[2];
    const Scalar initialAmplitude = amplitude;
    logFile << t0 << "," << amplitude << ",,," << rowDeviation << "\n";

    Scalar flux = 0.0, exactFlux = 0.0, maxRowDeviation = rowDeviation;
    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        const Scalar t = t0 + timeLoop->time()/timeScale;
        const Scalar tOld = t - timeLoop->previousTimeStepSize()/timeScale;
        const Scalar massOld = mass;
        measured = measure(x);
        amplitude = measured[0]; mass = measured[1]; rowDeviation = measured[2];
        // solute flux through the top from the change of the solute mass, in units of U
        flux = (mass - massOld)/(t - tOld);
        exactFlux = (exactMass(t) - exactMass(tOld))/(t - tOld);
        maxRowDeviation = std::max(maxRowDeviation, rowDeviation);
        logFile << t << "," << amplitude << "," << flux << "," << exactFlux << "," << rowDeviation << "\n";

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(dt);
    } while (!timeLoop->finished());

    vtkWriter.write(timeLoop->time());
    timeLoop->finalize(leafGridView.comm());

    std::ofstream profileFile(problem->name() + "_profile.csv");
    profileFile << std::setprecision(12) << "depth,amplitude\n";
    for (const auto& [z, row] : rowsOf(x))
        profileFile << (gridGeometry->bBoxMax()[1] - z)/height << "," << rowAmplitude(row) << "\n";

    const Scalar fluxError = std::abs(flux/exactFlux - 1.0);
    const Scalar logAmplification = seeded ? std::log(amplitude/initialAmplitude) : 0.0;
    std::cout << std::setprecision(10)
              << "Ra = " << rayleigh << ", k = " << k
              << ", log amplification " << logAmplification
              << ", relative flux error " << fluxError
              << ", max row deviation " << maxRowDeviation << std::endl;

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    if (hasParam("Problem.ReferenceLogAmplification"))
    {
        const Scalar error = std::abs(logAmplification - getParam<Scalar>("Problem.ReferenceLogAmplification"));
        std::cout << "error of the log amplification " << error << std::endl;
        if (!(error < getParam<Scalar>("Problem.LogAmplificationTolerance")))
            DUNE_THROW(Dune::Exception, "Error " << error << " of the log amplification exceeds the tolerance");
    }

    if (hasParam("Problem.FluxTolerance") && !(fluxError < getParam<Scalar>("Problem.FluxTolerance")))
        DUNE_THROW(Dune::Exception, "Relative flux error " << fluxError << " exceeds the tolerance");

    if (hasParam("Problem.MaxRowDeviationTolerance") && !(maxRowDeviation < getParam<Scalar>("Problem.MaxRowDeviationTolerance")))
        DUNE_THROW(Dune::Exception, "Horizontal variation " << maxRowDeviation << " exceeds the tolerance");

    return 0;
}
