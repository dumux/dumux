// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Onset of one-sided convection: growth of a perturbation of the diffusive boundary layer
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <vector>

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

    // amplitude of cos(kx) on each row of vertices of the tensor-product grid
    struct Row
    {
        Scalar volume = 0.0, mass = 0.0, cosine = 0.0, massCosine = 0.0, cosineSquared = 0.0;
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
            }
        }
        return rows;
    };

    const auto rowAmplitude = [](const Row& row)
    { return (row.massCosine - row.mass/row.volume*row.cosine)/row.cosineSquared; };

    // L2 norm of the mode amplitude over the depth
    const auto norm = [&](const SolutionVector& sol)
    {
        Scalar normSquared = 0.0;
        for (const auto& [z, row] : rowsOf(sol))
            normSquared += rowAmplitude(row)*rowAmplitude(row)*row.volume/(length*height);
        return std::sqrt(normSquared);
    };

    // maximum over the depth of the mode amplitude, as used for the growth rate by Riaz et al. (2006): the parabola through
    // the row of the maximum and its two neighbors gives the maximum between the rows
    const auto peak = [&](const SolutionVector& sol)
    {
        std::vector<Scalar> zRows, amplitudes;
        for (const auto& [zRow, row] : rowsOf(sol))
        {
            zRows.push_back(zRow);
            amplitudes.push_back(std::abs(rowAmplitude(row)));
        }
        const auto i = std::distance(amplitudes.begin(), std::max_element(amplitudes.begin(), amplitudes.end()));
        if (i == 0 || i == static_cast<std::ptrdiff_t>(amplitudes.size()) - 1)
            return amplitudes[i];
        const Scalar slopeLeft = (amplitudes[i] - amplitudes[i-1])/(zRows[i] - zRows[i-1]);
        const Scalar slopeRight = (amplitudes[i+1] - amplitudes[i])/(zRows[i+1] - zRows[i]);
        const Scalar curvature = (slopeRight - slopeLeft)/(zRows[i+1] - zRows[i-1]);
        if (curvature >= 0.0)
            return amplitudes[i];
        const Scalar zMax = 0.5*(zRows[i-1] + zRows[i]) - slopeLeft/(2.0*curvature);
        if (zMax < zRows[i-1] || zMax > zRows[i+1])
            return amplitudes[i];
        return amplitudes[i-1] + slopeLeft*(zMax - zRows[i-1]) + curvature*(zMax - zRows[i-1])*(zMax - zRows[i]);
    };

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
    logFile << std::setprecision(12) << "time,amplitude,maxAmplitude\n";
    Scalar amplitude = norm(x);
    const Scalar initialAmplitude = amplitude;
    logFile << t0 << "," << amplitude << "," << peak(x) << "\n";

    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        const Scalar t = t0 + timeLoop->time()/timeScale;
        amplitude = norm(x);
        logFile << t << "," << amplitude << "," << peak(x) << "\n";

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(dt);
    } while (!timeLoop->finished());

    vtkWriter.write(timeLoop->time());
    timeLoop->finalize(leafGridView.comm());

    std::ofstream profileFile(problem->name() + "_profile.csv");
    profileFile << std::setprecision(12) << "depth,amplitude\n";
    for (const auto& [z, row] : rowsOf(x))
        profileFile << (gridGeometry->bBoxMax()[1] - z)/height << "," << rowAmplitude(row) << "\n";

    const Scalar logAmplification = std::log(amplitude/initialAmplitude);
    std::cout << std::setprecision(10)
              << "Ra = " << rayleigh << ", k = " << k
              << ", log amplification " << logAmplification << std::endl;

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    if (hasParam("Problem.ReferenceLogAmplification"))
    {
        const Scalar error = std::abs(logAmplification - getParam<Scalar>("Problem.ReferenceLogAmplification"));
        std::cout << "error of the log amplification " << error << std::endl;
        if (!(error < getParam<Scalar>("Problem.LogAmplificationTolerance")))
            DUNE_THROW(Dune::Exception, "Error " << error << " of the log amplification exceeds the tolerance");
    }

    return 0;
}
