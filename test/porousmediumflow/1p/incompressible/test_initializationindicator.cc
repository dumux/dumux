// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Test the grid adaptation initialization indicator for the discretization given by TYPETAG
 *
 * The problem has Dirichlet boundaries at the bottom and the top, zero-flux Neumann
 * boundaries on the sides, and no sources. Exactly the elements with an intersection
 * on the Dirichlet boundaries have to be marked.
 */
#include <config.h>

#include <iostream>
#include <memory>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/adaptive/initializationindicator.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    initialize(argc, argv);
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

    GridAdaptInitializationIndicator<TypeTag> indicator(problem, gridGeometry, gridVariables);
    indicator.calculate(x);

    constexpr double eps = 1e-6;
    const auto yMin = gridGeometry->bBoxMin()[1];
    const auto yMax = gridGeometry->bBoxMax()[1];

    int numWrong = 0;
    int numMarked = 0;
    for (const auto& element : elements(leafGridView))
    {
        bool onDirichletBoundary = false;
        for (const auto& intersection : intersections(leafGridView, element))
        {
            const auto y = intersection.geometry().center()[1];
            if (intersection.boundary() && (y < yMin + eps || y > yMax - eps))
                onDirichletBoundary = true;
        }

        const int mark = indicator(element);
        numMarked += mark;
        if (mark != static_cast<int>(onDirichletBoundary))
        {
            std::cerr << "Element at " << element.geometry().center()
                      << (mark ? " is marked but has no Dirichlet boundary" : " has a Dirichlet boundary but is not marked")
                      << std::endl;
            ++numWrong;
        }
    }

    std::cout << "Marked " << numMarked << " elements, " << numWrong << " wrong" << std::endl;
    return numWrong == 0 ? 0 : 1;
}
