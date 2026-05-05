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

    using RotationTypeTag = Properties::TTag::MRPlateTestRotation;
    using DeformationTypeTag = Properties::TTag::MRPlateTestDeformation;
    using CommonTypeTag = Properties::TTag::MRPlateTestCommon;

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

    // Compute max element diameter (mesh size) and relative L2-error in vertical deformation w
    double hMax = 0.0;
    double l2norm2(0.0), l2norm2Ref(0.0);
    const auto& gg = *deformationGridGeometry;
    for (const auto& element : elements(gg.gridView()))
    {
        const auto geometry = element.geometry();
        hMax = std::max(hMax, Dumux::diameter(geometry));
        const auto elemSol = elementSolution(element, x[deformationIdx], gg);
        const auto& quad = Dune::QuadratureRules<double, Grid::dimension>::rule(geometry.type(), 3);
        for (const auto& qp : quad)
        {
            const auto globalPos = geometry.global(qp.position());
            const auto w = evalSolution(element, geometry, gg, elemSol, globalPos)[1];
            const auto wExact = deformationProblem->analyticDeformation(globalPos);
            const auto error = w - wExact;
            l2norm2 += error*error*qp.weight()*geometry.integrationElement(qp.position());
            l2norm2Ref += wExact*wExact*qp.weight()*geometry.integrationElement(qp.position());
        }
    }

    const auto relL2Error = std::sqrt(l2norm2/l2norm2Ref);
    std::cout << "Max element diameter: " << hMax << std::endl;
    std::cout << "Relative L2-error deformation w: " << relL2Error << std::endl;

    return 0;
}
