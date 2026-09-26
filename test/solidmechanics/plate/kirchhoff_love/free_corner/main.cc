// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/geometry/diameter.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using RotationTypeTag = Properties::TTag::KLFreeCornerTestRotation;
    using DeformationTypeTag = Properties::TTag::KLFreeCornerTestDeformation;
    using CommonTypeTag = Properties::TTag::KLFreeCornerTestCommon;

    using Grid = GetPropType<CommonTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using RotationGridGeometry = GetPropType<RotationTypeTag, Properties::GridGeometry>;
    using DeformationGridGeometry = GetPropType<DeformationTypeTag, Properties::GridGeometry>;
    auto rotationGridGeometry = std::make_shared<RotationGridGeometry>(leafGridView);
    auto deformationGridGeometry = std::make_shared<DeformationGridGeometry>(leafGridView);

    using CouplingManager = GetPropType<RotationTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(rotationGridGeometry, deformationGridGeometry);

    using RotationProblem = GetPropType<RotationTypeTag, Properties::Problem>;
    auto rotationProblem = std::make_shared<RotationProblem>(rotationGridGeometry, couplingManager);

    using DeformationProblem = GetPropType<DeformationTypeTag, Properties::Problem>;
    auto deformationProblem = std::make_shared<DeformationProblem>(deformationGridGeometry, couplingManager);

    constexpr auto rotationIdx = CouplingManager::rotationIdx;
    constexpr auto deformationIdx = CouplingManager::deformationIdx;
    using Traits = MultiDomainTraits<RotationTypeTag, DeformationTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    rotationProblem->applyInitialSolution(x[rotationIdx]);
    deformationProblem->applyInitialSolution(x[deformationIdx]);

    using RotationGridVariables = GetPropType<RotationTypeTag, Properties::GridVariables>;
    auto rotationGridVariables = std::make_shared<RotationGridVariables>(rotationProblem, rotationGridGeometry);

    using DeformationGridVariables = GetPropType<DeformationTypeTag, Properties::GridVariables>;
    auto deformationGridVariables = std::make_shared<DeformationGridVariables>(deformationProblem, deformationGridGeometry);

    couplingManager->init(rotationProblem, deformationProblem, x);
    rotationGridVariables->init(x[rotationIdx]);
    deformationGridVariables->init(x[deformationIdx]);

    VtkOutputModule<DeformationGridVariables, std::tuple_element_t<1, SolutionVector>>
        vtkWriter(*deformationGridVariables, x[deformationIdx], deformationProblem->name());
    vtkWriter.addVolumeVariable([](const auto& v){ return v.verticalDeformation(); }, "w");
    vtkWriter.addVolumeVariable([](const auto& v){ return v.shearCurlPotential(); }, "psi");
    vtkWriter.addVolumeVariable([](const auto& v){ return v.shearGradPotential(); }, "phi");
    vtkWriter.write(0.0);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(rotationProblem, deformationProblem),
        std::make_tuple(rotationGridGeometry, deformationGridGeometry),
        std::make_tuple(rotationGridVariables, deformationGridVariables),
        couplingManager
    );

    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);
    nonLinearSolver->solve(x);

    vtkWriter.write(1.0);

    double hMax = 0.0;
    const auto& gg = *deformationGridGeometry;
    for (const auto& element : elements(gg.gridView()))
        hMax = std::max(hMax, Dumux::diameter(element.geometry()));

    const auto nu = getParam<double>("Problem.PoissonRatio");
    const auto E = getParam<double>("Problem.E");
    const auto t = getParam<double>("Problem.Thickness");
    const auto D = E*t*t*t/(12.0*(1.0 - nu*nu));
    const auto length = getParam<double>("Problem.Length", 1.0);
    using Position = Dune::FieldVector<double, 2>;

    const auto normalsAt = [&](const Position& point)
    {
        std::vector<Position> normals;
        for (const auto& element : elements(gg.gridView()))
            for (const auto& is : intersections(gg.gridView(), element))
            {
                if (!is.boundary())
                    continue;
                const auto geo = is.geometry();
                bool touches = false;
                for (int c = 0; c < geo.corners(); ++c)
                    if ((geo.corner(c) - point).two_norm() < 1e-8*length)
                        touches = true;
                if (!touches)
                    continue;
                const auto n = is.centerUnitOuterNormal();
                if (std::none_of(normals.begin(), normals.end(),
                                 [&](const auto& m){ return (m - n).two_norm() < 1e-6; }))
                    normals.push_back(n);
            }
        return normals;
    };

    // Element-wise rotation gradients need not agree at a shared vertex.
    const auto momentAt = [&](const Position& point)
    {
        Dune::FieldMatrix<double, 2, 2> M(0.0);
        int count = 0;
        for (const auto& element : elements(rotationGridGeometry->gridView()))
        {
            const auto geo = element.geometry();
            bool touches = false;
            for (int c = 0; c < geo.corners(); ++c)
                if ((geo.corner(c) - point).two_norm() < 1e-8*length)
                    touches = true;
            if (!touches)
                continue;
            const auto elemSol = elementSolution(element, x[rotationIdx], *rotationGridGeometry);
            const auto g = evalGradients(element, geo, *rotationGridGeometry, elemSol, point);
            M[0][0] += -D*(g[0][0] + nu*g[1][1]);
            M[1][1] += -D*(g[1][1] + nu*g[0][0]);
            M[0][1] += -D*(1.0 - nu)*0.5*(g[0][1] + g[1][0]);
            ++count;
        }
        if (count > 0)
            M /= count;
        M[1][0] = M[0][1];
        return M;
    };

    const auto cornerForce = [&](const Position& point)
    {
        const auto normals = normalsAt(point);
        if (normals.size() != 2)
            return std::numeric_limits<double>::quiet_NaN();
        const auto M = momentAt(point);
        const auto twisting = [&](Position n)
        {
            Position s{n[1], -n[0]}; // s = J n
            Position Mn(0.0);
            M.mv(n, Mn);
            return s*Mn;
        };
        return twisting(normals[0]) - twisting(normals[1]);
    };

    double wMax = 0.0;
    for (const auto& dof : x[deformationIdx])
        wMax = std::max(wMax, std::abs(dof[1]));
    const auto scale = D*(1.0 - nu)*wMax/(length*length);

    std::cout << "Max element diameter: " << hMax << std::endl;
    std::cout << "Maximum deflection: " << wMax << std::endl;

    const auto corners = getParam<std::vector<double>>("Problem.Corners", std::vector<double>{});
    for (std::size_t i = 0; 2*i + 1 < corners.size(); ++i)
    {
        Position point{corners[2*i], corners[2*i + 1]};
        const auto clamped = deformationProblem->onClampedEdge(point);
        const auto force = cornerForce(point);
        std::cout << "Corner (" << point[0] << "," << point[1] << ") "
                  << (clamped ? "clamped-free" : "free-free   ")
                  << "  R = " << force
                  << "  ratio = " << std::abs(force)/scale << std::endl;
    }

    return 0;
}
