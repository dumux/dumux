// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Test that control-volume schemes are exact for flow through two layers of different permeability
 *
 * The interface between the layers coincides with element faces. With a constant permeability per element,
 * the Darcy velocity of the exact pressure is constant and identical in all elements, so the exact pressure
 * makes the residual of every control volume vanish.
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
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/discretization/cvfe/interpolate.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

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

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    std::vector<double> coefficients(gridGeometry->numDofs());
    CVFE::interpolate(*gridGeometry, coefficients, [&](const auto& globalPos) { return problem->exactPressure(globalPos); });
    for (std::size_t dofIdx = 0; dofIdx < x.size(); ++dofIdx)
        x[dofIdx] = coefficients[dofIdx];

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

    // the residual of a control volume sums the fluxes through faces of size h^(dim-1),
    // and the Darcy velocity is one
    using GlobalPosition = typename Grid::template Codim<0>::Geometry::GlobalCoordinate;
    constexpr int dim = Grid::dimension;
    const auto upperRight = getParam<GlobalPosition>("Grid.UpperRight");
    const auto cells = getParam<std::array<int, dim>>("Grid.Cells");
    const double faceSize = std::pow(upperRight[0]/cells[0], dim - 1);

    double maxResidual = 0.0;
    std::size_t numChecked = 0;
    for (std::size_t dofIdx = 0; dofIdx < x.size(); ++dofIdx)
    {
        if (isConstrained[dofIdx])
            continue;
        maxResidual = std::max(maxResidual, std::abs(assembler->residual()[dofIdx][0])/faceSize);
        ++numChecked;
    }

    std::cout << "max |residual|/h^(dim-1) over " << numChecked << " unconstrained dofs: " << maxResidual << std::endl;

    if (numChecked == 0)
        DUNE_THROW(Dune::Exception, "No unconstrained dofs");
    if (!(maxResidual < 1e-10))
        DUNE_THROW(Dune::Exception, "The exact pressure does not satisfy the discrete equations");

    return 0;
}
