// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Hydrostatic test for the one-phase model on a grid whose cell-centre connections
 *        are not parallel to the face normals. The discrete solution has to reproduce the
 *        hydrostatic pressure and vanishing fluxes (up to round-off).
 */

#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/pdesolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager_ug.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 1.0e5;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    LinearPDESolver<Assembler, LinearSolver> solver(assembler, linearSolver);
    solver.solve(x);
    gridVariables->update(x);

    // scales for relative errors: hydrostatic pressure difference over the domain height
    // and the Darcy velocity magnitude of the gravity-driven flux
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    const auto& bBoxMin = gridGeometry->bBoxMin();
    const auto& bBoxMax = gridGeometry->bBoxMax();
    const Scalar T = problem->spatialParams().temperatureAtPos(bBoxMin);
    const Scalar rho = FluidSystem::density(T, 1.0e5);
    const Scalar mu = FluidSystem::viscosity(T, 1.0e5);
    const Scalar g = problem->spatialParams().gravity(bBoxMin).two_norm();
    const Scalar pressureScale = rho*g*(bBoxMax[1] - bBoxMin[1]);
    const Scalar kMax = std::max(getParam<Scalar>("SpatialParams.PermeabilityLeft"),
                                 getParam<Scalar>("SpatialParams.PermeabilityRight"));
    const Scalar velocityScale = kMax/mu*rho*g;

    Scalar maxPressureError = 0.0;
    Scalar maxFluxError = 0.0;
    auto fvGeometry = localView(*gridGeometry);
    auto elemVolVars = localView(gridVariables->curGridVolVars());
    auto elemFluxVarsCache = localView(gridVariables->gridFluxVarsCache());
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    const auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
    for (const auto& element : elements(leafGridView))
    {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, x);
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        for (const auto& scv : scvs(fvGeometry))
        {
            const Scalar error = std::abs(x[scv.dofIndex()][0] - problem->exactPressure(scv.dofPosition()));
            std::cout << "dof at (" << scv.dofPosition() << "): p = " << x[scv.dofIndex()][0]
                      << ", exact " << problem->exactPressure(scv.dofPosition()) << ", error " << error << " Pa\n";
            maxPressureError = std::max(maxPressureError, error/pressureScale);
        }

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
                continue;
            FluxVariables fluxVars;
            fluxVars.init(*problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
            const Scalar flux = fluxVars.advectiveFlux(0, upwindTerm);
            maxFluxError = std::max(maxFluxError, std::abs(flux)/(velocityScale*scvf.area()));
        }
    }

    std::cout << "maximum relative pressure error: " << maxPressureError << "\n"
              << "maximum relative interior flux: " << maxFluxError << std::endl;

    if (maxPressureError > 1e-10 || maxFluxError > 1e-10)
    {
        std::cerr << "The hydrostatic solution is not reproduced." << std::endl;
        return 1;
    }

    if (Dune::MPIHelper::getCommunication().rank() == 0)
        Parameters::print();

    return 0;
}
