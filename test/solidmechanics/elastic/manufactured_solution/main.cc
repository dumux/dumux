// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Test for the linear elastic model with a manufactured solution
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#if DUMUX_HAVE_GRIDFORMAT
#include <dumux/io/gridwriter.hh>
#include <dumux/io/cvfegridfunction.hh>
#endif

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TestElastic;

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

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    // L2 error of the displacement
    double l2Error = 0.0;
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(leafGridView))
    {
        elemDisc.bind(element);
        const auto geometry = element.geometry();
        const auto& localBasis = elemDisc.feLocalBasis();
        using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
        for (const auto& qp : Dune::QuadratureRules<double, 2>::rule(geometry.type(), 6))
        {
            const auto global = geometry.global(qp.position());
            const FEInterpolationPointData<GlobalPosition, std::decay_t<decltype(localBasis)>>
                shapeData(geometry, qp.position(), global, localBasis);
            Dune::FieldVector<double, 2> u(0.0);
            for (const auto& localDof : localDofs(elemDisc))
                u.axpy(shapeData.shapeValues()[localDof.index()][0], x[localDof.dofIndex()]);
            u -= problem->exactSolution(global);
            l2Error += u.two_norm2()*qp.weight()*geometry.integrationElement(qp.position());
        }
    }
    std::cout << "L2 error displacement: " << std::sqrt(l2Error) << std::endl;

#if DUMUX_HAVE_GRIDFORMAT
    SolutionVector xExact(x.size());
    for (const auto& element : elements(leafGridView))
    {
        elemDisc.bind(element);
        for (const auto& localDof : localDofs(elemDisc))
            xExact[localDof.dofIndex()] = problem->exactSolution(ipData(elemDisc, localDof).global());
    }

    IO::GridWriter writer{IO::Format::vtu, leafGridView, IO::order<1>};
    writer.setPointField("u", IO::cvfeGridFunction(*gridGeometry, x));
    writer.setPointField("u_exact", IO::cvfeGridFunction(*gridGeometry, xExact));
    writer.write(problem->name());
#endif

    return 0;
}
