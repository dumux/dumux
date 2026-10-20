// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Compares the Jacobian of the second-order shallow water model assembled with the two-ring stencil
 *        with a central finite-difference Jacobian of the global residual, at a smooth wet state on the
 *        lake-at-rest geometry, for the stationary assembler and for the multi-stage assembler of the time
 *        integration methods (residual of an implicit solve with storage and spatial operator). A dependency
 *        missing from the connectivity map, or a race between elements assembled concurrently, shows up as a
 *        difference of the order of the entries. With the limiter, its kinks cause small differences, so the
 *        tolerance is a parameter.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/experimental/assembly/multistagefvassembler.hh>
#include <dumux/experimental/timestepping/multistepmethods.hh>
#include <dumux/experimental/timestepping/timelevel.hh>

#include "../lakeatrest/properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::LakeAtRestSecondOrder;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // smooth wet state that is not at rest
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    for (const auto& element : elements(leafGridView))
    {
        const auto eIdx = gridGeometry->elementMapper().index(element);
        const auto pos = element.geometry().center();
        const Scalar bed = problem->spatialParams().bedSurfaceAtPos(pos);
        x[eIdx][0] = 1.0 - bed + 0.05*std::sin(2.0*pos[0])*std::cos(3.0*pos[1]);
        x[eIdx][1] = 0.3*std::cos(1.5*pos[0] + pos[1]);
        x[eIdx][2] = -0.2*std::sin(2.5*pos[1] - 0.5*pos[0]);
    }

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // face neighbors, to show that the two-ring entries are needed
    std::vector<std::set<std::size_t>> neighbors(gridGeometry->numDofs());
    for (const auto& element : elements(leafGridView))
        for (const auto& intersection : intersections(leafGridView, element))
            if (intersection.neighbor())
                neighbors[gridGeometry->elementMapper().index(element)].insert(gridGeometry->elementMapper().index(intersection.outside()));

    const auto tolerance = getParam<Scalar>("Problem.JacobianTolerance");

    // compares an assembled Jacobian with central differences of the residual function at x
    auto checkJacobian = [&](const std::string& name, const auto& A, auto&& residual)
    {
        Scalar maxAbs = 0.0;
        for (auto row = A.begin(); row != A.end(); ++row)
            for (auto col = row->begin(); col != row->end(); ++col)
                maxAbs = std::max(maxAbs, col->infinity_norm());

        const Scalar delta = 1e-6;
        Scalar maxDifference = 0.0;
        std::size_t numBeyondFaceNeighbors = 0;
        for (std::size_t k = 0; k < x.size(); ++k)
        {
            for (int pvIdx = 0; pvIdx < 3; ++pvIdx)
            {
                auto evalResidual = [&](Scalar shift) {
                    auto xs = x;
                    xs[k][pvIdx] += shift;
                    gridVariables->update(xs);
                    return residual(xs);
                };

                const auto rPlus = evalResidual(delta);
                const auto rMinus = evalResidual(-delta);
                for (std::size_t i = 0; i < x.size(); ++i)
                {
                    for (int eqIdx = 0; eqIdx < 3; ++eqIdx)
                    {
                        const Scalar fd = (rPlus[i][eqIdx] - rMinus[i][eqIdx])/(2.0*delta);
                        const Scalar assembled = A.exists(i, k) ? A[i][k][eqIdx][pvIdx] : 0.0;
                        maxDifference = std::max(maxDifference, std::abs(fd - assembled));
                        if (i != k && !neighbors[i].count(k) && std::abs(fd) > 1e-6*maxAbs)
                            ++numBeyondFaceNeighbors;
                    }
                }
            }
        }
        gridVariables->update(x);

        std::cout << name << ": max |J_assembled - J_fd| / max |J| = " << maxDifference/maxAbs
                  << ", entries beyond face neighbors: " << numBeyondFaceNeighbors << std::endl;

        if (!(maxDifference < tolerance*maxAbs))
            DUNE_THROW(Dune::Exception, name << ": assembled Jacobian differs from the finite-difference Jacobian");
        if (numBeyondFaceNeighbors == 0)
            DUNE_THROW(Dune::Exception, name << ": no dependency beyond the face neighbors, the test does not check the two-ring stencil");
    };

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    Assembler assembler(problem, gridGeometry, gridVariables);
    assembler.assembleJacobianAndResidual(x);
    checkJacobian("Stationary assembler", assembler.jacobian(), [&](const SolutionVector& xs) {
        typename Assembler::ResidualType r(xs.size());
        r = 0.0;
        assembler.assembleResidual(r, xs);
        return r;
    });

    // residual w_M M + w_R R of an implicit solve, as assembled for the time integration methods
    using MultiStageAssembler = Experimental::MultiStageFVAssembler<TypeTag, DiffMethod::numeric>;
    using SolveParams = typename MultiStageAssembler::SolveParams;
    auto method = std::make_shared<Experimental::MultiStep::BDF<Scalar>>(2);
    const auto xOld = x;
    MultiStageAssembler multiStageAssembler(problem, gridGeometry, gridVariables, method, xOld);
    multiStageAssembler.setLinearSystem();
    const Experimental::TimeLevel<Scalar> timeLevel(0.0);
    const Scalar storageWeight = 1.0/getParam<Scalar>("Problem.TimeStepSize", 0.1);
    const Scalar spatialWeight = 1.0;
    multiStageAssembler.setImplicitSolve(x, SolveParams{timeLevel, {0.0, storageWeight, spatialWeight}, nullptr});
    multiStageAssembler.assembleJacobianAndResidual(x);
    checkJacobian("Multi-stage assembler", multiStageAssembler.jacobian(), [&](const SolutionVector& xs) {
        typename MultiStageAssembler::ResidualType storage, spatial;
        multiStageAssembler.assembleOperators(xs, timeLevel, {.storage = &storage, .spatial = &spatial});
        storage *= storageWeight;
        storage.axpy(spatialWeight, spatial);
        return storage;
    });

    return 0;
}
