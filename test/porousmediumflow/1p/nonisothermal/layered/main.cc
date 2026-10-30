// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Test that the box scheme with the face average of the thermal conductivity is exact for heat
 *        conduction through two layers
 *
 * The interface between the layers coincides with faces of the box control volumes, so the conductivity is
 * constant in each control volume. Between two neighboring control volumes, the heat flows through two
 * half-cells in series, whose conductivity is the harmonic mean of the two conductivities, as given by
 * `Flux.ScvfThermalConductivity = FaceTensorAverage`. The exact temperature therefore makes the residual
 * of every control volume vanish.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/discretization/cvfe/interpolate.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::OnePNILayeredBox;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    std::vector<double> temperature(gridGeometry->numDofs());
    CVFE::interpolate(*gridGeometry, temperature, [&](const auto& globalPos) { return problem->exactTemperature(globalPos); });
    for (std::size_t dofIdx = 0; dofIdx < x.size(); ++dofIdx)
    {
        x[dofIdx][Indices::pressureIdx] = problem->pressure();
        x[dofIdx][Indices::temperatureIdx] = temperature[dofIdx];
    }

    std::vector<bool> isConstrained(gridGeometry->numDofs(), false);
    for (const auto& constraint : problem->constraints())
        isConstrained[constraint.dofIndex()] = true;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    assembler->setLinearSystem();
    assembler->assembleResidual(x);

    // the residual of a control volume sums the heat fluxes through faces of size h^(dim-1)
    using GlobalPosition = typename Grid::template Codim<0>::Geometry::GlobalCoordinate;
    constexpr int dim = Grid::dimension;
    const auto upperRight = getParam<GlobalPosition>("Grid.UpperRight");
    const auto cells = getParam<std::array<int, dim>>("Grid.Cells");
    const double faceSize = std::pow(upperRight[0]/cells[0], dim - 1);
    const double scale = std::abs(problem->heatFlux())*faceSize;

    double maxResidual = 0.0;
    std::size_t numChecked = 0;
    for (std::size_t dofIdx = 0; dofIdx < x.size(); ++dofIdx)
    {
        if (isConstrained[dofIdx])
            continue;
        maxResidual = std::max(maxResidual, std::abs(assembler->residual()[dofIdx][Indices::energyEqIdx])/scale);
        ++numChecked;
    }

    std::cout << "max |energy residual|/(|q| h^(dim-1)) over " << numChecked << " unconstrained dofs: " << maxResidual << std::endl;

    if (numChecked == 0)
        DUNE_THROW(Dune::Exception, "No unconstrained dofs");
    if (!(maxResidual < 1e-10))
        DUNE_THROW(Dune::Exception, "The exact temperature does not satisfy the discrete equations");

    return 0;
}
