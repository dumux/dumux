// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Growing layered shell: the growth is increased in steps and the L2 and H1 errors of the
 *        displacement against the radially symmetric solution are reported at the final growth
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"
#include "shellgrid.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    static constexpr int dim = Grid::dimension;
    const std::array<Scalar, 3> radii{getParam<Scalar>("Problem.InnerRadius"),
                                      getParam<Scalar>("Problem.InterfaceRadius"),
                                      getParam<Scalar>("Problem.OuterRadius")};
    const auto cellParam = getParam<std::vector<int>>("Grid.Cells");
    const auto grid = createLayeredShellGrid<Grid>(radii, std::array<int, 3>{cellParam[0], cellParam[1], cellParam[2]});
    const auto& leafGridView = grid->leafGridView();

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

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonlinearSolver(assembler, linearSolver);

    // the growth is increased in steps, each starting from the previous solution
    const auto numSteps = getParam<int>("Problem.GrowthSteps");
    for (int step = 1; step <= numSteps; ++step)
    {
        problem->setGrowthFraction(static_cast<Scalar>(step)/numSteps);
        nonlinearSolver.solve(x);
    }

    if (getParam<bool>("Problem.EnableVtkOutput", false))
    {
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::box)
        {
            VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
            vtkWriter.addField(x, "d");
            vtkWriter.write(1.0);
        }
    }

    Scalar l2Error = 0.0, h1Error = 0.0, hMax = 0.0;
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(leafGridView))
    {
        elemDisc.bind(element);
        const auto& geometry = element.geometry();
        for (int i = 0; i < geometry.corners(); ++i)
            for (int j = i + 1; j < geometry.corners(); ++j)
                hMax = std::max(hMax, (geometry.corner(i) - geometry.corner(j)).two_norm());

        const auto& localBasis = elemDisc.feLocalBasis();
        using ShapeData = FEInterpolationPointData<typename Grid::template Codim<0>::Geometry::GlobalCoordinate,
                                                   std::decay_t<decltype(localBasis)>>;
        for (const auto& qp : Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), 8))
        {
            const auto global = geometry.global(qp.position());
            const ShapeData shapeData(geometry, qp.position(), global, localBasis);
            Dune::FieldVector<Scalar, dim> u(0.0);
            Dune::FieldMatrix<Scalar, dim, dim> gradU(0.0);
            for (const auto& localDof : localDofs(elemDisc))
            {
                const auto& dofValue = x[localDof.dofIndex()];
                u.axpy(shapeData.shapeValues()[localDof.index()][0], dofValue);
                for (int d = 0; d < dim; ++d)
                    gradU[d].axpy(dofValue[d], shapeData.gradN(localDof.index()));
            }
            u -= problem->exactDisplacement(global);
            gradU -= problem->exactDisplacementGradient(global);
            const Scalar weight = qp.weight()*geometry.integrationElement(qp.position());
            l2Error += weight*u.two_norm2();
            h1Error += weight*gradU.frobenius_norm2();
        }
    }

    std::cout << "Max element diameter: " << hMax << "\n"
              << "L2 error displacement: " << std::sqrt(l2Error) << "\n"
              << "H1 error displacement: " << std::sqrt(h1Error) << "\n";

    return 0;
}
