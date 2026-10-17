// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Growth bending strip in the mixed displacement-pressure formulation: the curvature is
 *        increased in steps and the errors of displacement and pressure against
 *        the exact arc are reported at the final curvature
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/traits.hh>

#include "properties_mixed.hh"
#include <dumux/io/grid/gridmanager_yasp.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;
    using MomTypeTag = Properties::TTag::MOMENTUM_TYPETAG;
    using PresTypeTag = Properties::TTag::BendingStripPressure;
    using MDTraits = MultiDomainTraits<MomTypeTag, PresTypeTag>;
    using Scalar = typename MDTraits::Scalar;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<MomTypeTag, Properties::Grid>;
    static constexpr int dim = Grid::dimension;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using MomGG = GetPropType<MomTypeTag, Properties::GridGeometry>;
    using PresGG = GetPropType<PresTypeTag, Properties::GridGeometry>;
    auto momGG = std::make_shared<MomGG>(leafGridView);
    auto presGG = std::make_shared<PresGG>(leafGridView);

    using CouplingManager = GetPropType<MomTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(momGG, presGG);

    using MomProblem = GetPropType<MomTypeTag, Properties::Problem>;
    using PresProblem = GetPropType<PresTypeTag, Properties::Problem>;
    auto momProblem = std::make_shared<MomProblem>(momGG, couplingManager);
    auto presProblem = std::make_shared<PresProblem>(presGG, couplingManager);

    using SolutionVector = typename MDTraits::SolutionVector;
    SolutionVector x;
    x[CouplingManager::momentumIdx].resize(momGG->numDofs());
    x[CouplingManager::pressureIdx].resize(presGG->numDofs());
    x = 0.0;

    using MomGV = GetPropType<MomTypeTag, Properties::GridVariables>;
    using PresGV = GetPropType<PresTypeTag, Properties::GridVariables>;
    auto momGV = std::make_shared<MomGV>(momProblem, momGG);
    auto presGV = std::make_shared<PresGV>(presProblem, presGG);
    momGV->init(x[CouplingManager::momentumIdx]);
    presGV->init(x[CouplingManager::pressureIdx]);

    couplingManager->init(momProblem, presProblem, x);
    couplingManager->computeColorsForAssembly();

    using Assembler = Experimental::MultiDomainAssembler<MDTraits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momProblem, presProblem),
                                                 std::make_tuple(momGG, presGG),
                                                 std::make_tuple(momGV, presGV),
                                                 couplingManager);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager> nonlinearSolver(assembler, linearSolver, couplingManager);

    // the curvature is increased in steps, each starting from the previous solution
    const auto curvature = getParam<Scalar>("Problem.Curvature");
    const auto numSteps = getParam<int>("Problem.GrowthSteps");
    for (int step = 1; step <= numSteps; ++step)
    {
        momProblem->setCurvature(curvature*step/numSteps);
        presProblem->setCurvature(curvature*step/numSteps);
        nonlinearSolver.solve(x);
    }

    Scalar l2Error = 0.0, h1Error = 0.0, pressureError = 0.0, hMax = 0.0;
    auto momElemDisc = localView(*momGG);
    auto presElemDisc = localView(*presGG);
    const auto& u = x[CouplingManager::momentumIdx];
    const auto& p = x[CouplingManager::pressureIdx];
    for (const auto& element : elements(leafGridView))
    {
        momElemDisc.bind(element);
        presElemDisc.bind(element);
        const auto& geometry = element.geometry();
        for (int i = 0; i < geometry.corners(); ++i)
            for (int j = i + 1; j < geometry.corners(); ++j)
                hMax = std::max(hMax, (geometry.corner(i) - geometry.corner(j)).two_norm());

        const auto& momBasis = momElemDisc.feLocalBasis();
        const auto& presBasis = presElemDisc.feLocalBasis();
        using GlobalPosition = typename Grid::template Codim<0>::Geometry::GlobalCoordinate;
        using MomShapeData = FEInterpolationPointData<GlobalPosition, std::decay_t<decltype(momBasis)>>;
        using PresShapeData = FEInterpolationPointData<GlobalPosition, std::decay_t<decltype(presBasis)>>;
        for (const auto& qp : Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), 8))
        {
            const auto global = geometry.global(qp.position());
            const MomShapeData momShape(geometry, qp.position(), global, momBasis);
            const PresShapeData presShape(geometry, qp.position(), global, presBasis);

            Dune::FieldVector<Scalar, dim> uh(0.0);
            Dune::FieldMatrix<Scalar, dim, dim> gradUh(0.0);
            for (const auto& localDof : localDofs(momElemDisc))
            {
                const auto& dofValue = u[localDof.dofIndex()];
                uh.axpy(momShape.shapeValues()[localDof.index()][0], dofValue);
                for (int d = 0; d < dim; ++d)
                    gradUh[d].axpy(dofValue[d], momShape.gradN(localDof.index()));
            }
            Scalar ph = 0.0;
            for (const auto& localDof : localDofs(presElemDisc))
                ph += presShape.shapeValues()[localDof.index()][0]*p[localDof.dofIndex()][0];

            uh -= momProblem->growth().exactDisplacement(global);
            gradUh -= momProblem->growth().exactDisplacementGradient(global);
            const Scalar weight = qp.weight()*geometry.integrationElement(qp.position());
            l2Error += weight*uh.two_norm2();
            h1Error += weight*gradUh.frobenius_norm2();
            pressureError += weight*ph*ph;
        }
    }

    std::cout << "Max element diameter: " << hMax << "\n"
              << "L2 error displacement: " << std::sqrt(l2Error) << "\n"
              << "H1 error displacement: " << std::sqrt(h1Error) << "\n"
              << "L2 error pressure: " << std::sqrt(pressureError) << "\n";

    return 0;
}
