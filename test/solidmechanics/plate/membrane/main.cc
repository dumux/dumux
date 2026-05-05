// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <iostream>
#include <cmath>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/geometry/diameter.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::MembranePlateTest;
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
    problem->applyInitialSolution(x);

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addVolumeVariable([](const auto& v){ return v.deformation(); }, "w");
    vtkWriter.write(0.0);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    gridVariables->update(x);
    vtkWriter.write(1.0);

    // Compute max element diameter (mesh size) and relative L2-error in vertical deformation w
    double hMax = 0.0;
    double l2norm2(0.0), l2norm2Ref(0.0);
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        const auto geometry = element.geometry();
        hMax = std::max(hMax, Dumux::diameter(geometry));
        const auto elemSol = elementSolution(element, x, *gridGeometry);
        const auto& quad = Dune::QuadratureRules<double, Grid::dimension>::rule(geometry.type(), 3);
        for (const auto& qp : quad)
        {
            const auto globalPos = geometry.global(qp.position());
            const auto w = evalSolution(element, geometry, *gridGeometry, elemSol, globalPos)[0];
            const auto wExact = problem->analyticDeformation(globalPos);
            const auto error = w - wExact;
            l2norm2 += error*error*qp.weight()*geometry.integrationElement(qp.position());
            l2norm2Ref += wExact*wExact*qp.weight()*geometry.integrationElement(qp.position());
        }
    }

    std::cout << "Max element diameter: " << hMax << std::endl;
    std::cout << "Relative L2-error deformation w: " << std::sqrt(l2norm2/l2norm2Ref) << std::endl;

    return 0;
}
