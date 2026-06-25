// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Dimensionless Boussinesq dissolution test: one-sided Rayleigh-Bénard convection.
 */

#include <config.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    constexpr int dimW = GetPropType<TypeTag, Properties::Grid>::dimensionworld;
    const std::string defaultParamFile = (dimW == 3) ? "params3d.input" : "params.input";
    Parameters::init(argc, argv, [](auto&){}, defaultParamFile);
    GetPropType<TypeTag, Properties::FluidSystem>::init();

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);

    const Scalar Ra = getParam<Scalar>("DimensionlessNumbers.Ra");

    // Optionally start from the analytical diffusive profile at t = tStart.
    // Default: tStart = 146/Ra^2 (dimensionless critical time from linear stability).
    // Set TimeLoop.TStart = 0 in params.input to use the standard zero-concentration IC.
    const Scalar tcrit  = 146.0 / Ra;
    const Scalar tStart = getParam<Scalar>("TimeLoop.TStart", tcrit);
    if (tStart > 0.0)
    {
        const Scalar zMax = gridGeometry->bBoxMax()[dimWorld-1];
        for (const auto& vertex : vertices(leafGridView))
        {
            const auto idx = gridGeometry->dofMapper().index(vertex);
            const Scalar z   = vertex.geometry().center()[dimWorld-1];
            const Scalar eta = (zMax - z) / (2.0 * std::sqrt(tStart / Ra));
            x[idx][Indices::concentrationIdx] = std::erfc(eta);
        }
    }

    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // --- Sherwood number: dimensionless dissolution flux at the top boundary ---
    // Sh(t) = (1/A_top) * integral_top (dC/dn) dA  where n points into the domain
    // For pure diffusion: Sh = sqrt(Ra/(pi*t)); convection enhances this.
    const Scalar zMax = gridGeometry->bBoxMax()[dimWorld-1];

    // horizontal area of the top face (all dimensions except the vertical one)
    Scalar domainTopArea = 1.0;
    for (int d = 0; d < dimWorld-1; ++d)
        domainTopArea *= gridGeometry->bBoxMax()[d] - gridGeometry->bBoxMin()[d];

    auto computeSherwood = [&]() -> Scalar {
        Scalar totalFlux = 0.0;
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);

            bool hasTopFace = false;
            for (const auto& scvf : scvfs(fvGeometry))
                if (scvf.boundary() && scvf.center()[dimWorld-1] > zMax - 1e-8)
                    { hasTopFace = true; break; }
            if (!hasTopFace) continue;

            const auto elemSol = elementSolution(element, x, *gridGeometry);
            const auto geo = element.geometry();

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (!scvf.boundary() || scvf.center()[dimWorld-1] < zMax - 1e-8)
                    continue;

                const auto grads = evalGradients(element, geo, *gridGeometry, elemSol, scvf.center());
                totalFlux += (grads[Indices::concentrationIdx] * scvf.unitOuterNormal()) * scvf.area();
            }
        }
        // sum over all MPI ranks
        totalFlux = leafGridView.comm().sum(totalFlux);
        return totalFlux / domainTopArea;
    };

    // --- Mass-balance flux: integrate C over domain, use dM/dt as dissolution rate ---
    auto computeTotalMass = [&]() -> Scalar {
        Scalar mass = 0.0;
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                mass += x[scv.dofIndex()][Indices::concentrationIdx] * scv.volume();
        }
        return leafGridView.comm().sum(mass);
    };

    std::ofstream shFile(problem->name() + "_sherwood.csv");
    shFile << std::scientific << std::setprecision(6);
    // Sh_grad / F_grad: gradient at top interface (diffusive flux proxy)
    // Sh_mass / F_mass: dM/dt mass balance (total flux, diffusion + advection)
    shFile << "time,Sh_grad,F_grad,Sh_mass,F_mass\n";
    // -------------------------------------------------------------------------

    const auto tEnd       = getParam<Scalar>("TimeLoop.TEnd");
    const auto dtInit     = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt      = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto vtkInterval = getParam<Scalar>("Output.VtkOutputInterval", 10.0);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    const auto writeSherwood = [&](Scalar t, Scalar massFluxPerWidth = 0.0) {
        const Scalar Sh_grad = computeSherwood();
        const Scalar Sh_mass = massFluxPerWidth * Ra;
        shFile << t << "," << Sh_grad << "," << Sh_grad/Ra
                    << "," << Sh_mass << "," << massFluxPerWidth << "\n";
        shFile.flush();
    };

    Scalar massOld = computeTotalMass();
    writeSherwood(tStart);

    Scalar nextVtkTime = tStart + vtkInterval;
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(tStart, dtInit, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits,
                         LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();


    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);

        const Scalar dt = timeLoop->timeStepSize();
        const Scalar massNew = computeTotalMass();
        const Scalar massFluxPerWidth = (massNew - massOld) / (dt * domainTopArea);
        massOld = massNew;

        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        if (timeLoop->time() >= nextVtkTime - 1e-8*vtkInterval || timeLoop->finished())
        {
            vtkWriter.write(timeLoop->time());
            nextVtkTime += vtkInterval;
        }

        writeSherwood(timeLoop->time(), massFluxPerWidth);

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}