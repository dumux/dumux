// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup KirchhoffLovePlate
 * \brief Test that the Kirchhoff-Love plate coupling manager, used as a
 *        default-constructed sub-manager attached to an externally owned solution
 *        storage, assembles the same Jacobian and residual as the regular manager.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>

#include <dune/common/shared_ptr.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/io/grid/gridmanager_foam.hh>

#include "properties.hh"

namespace Dumux {

// deviation between two sparse matrices with (expected) identical sparsity pattern
template<class Matrix>
double maxAbsDifferenceMatrix(const Matrix& a, const Matrix& b)
{
    if (a.N() != b.N() || a.M() != b.M())
        return std::numeric_limits<double>::infinity();

    std::size_t nnzA = 0, nnzB = 0;
    double maxDiff = 0.0;
    for (auto rowIt = a.begin(); rowIt != a.end(); ++rowIt)
    {
        nnzA += rowIt->size();
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
        {
            if (!b.exists(rowIt.index(), colIt.index()))
                return std::numeric_limits<double>::infinity();

            const auto& blockB = b[rowIt.index()][colIt.index()];
            for (std::size_t i = 0; i < colIt->N(); ++i)
                for (std::size_t j = 0; j < colIt->M(); ++j)
                    maxDiff = std::max(maxDiff, std::abs((*colIt)[i][j] - blockB[i][j]));
        }
    }

    for (auto rowIt = b.begin(); rowIt != b.end(); ++rowIt)
        nnzB += rowIt->size();

    return nnzA == nnzB ? maxDiff : std::numeric_limits<double>::infinity();
}

// deviation between two block vectors of identical size
template<class BlockVector>
double maxAbsDifferenceVector(const BlockVector& a, const BlockVector& b)
{
    double maxDiff = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[i].size(); ++j)
            maxDiff = std::max(maxDiff, std::abs(a[i][j] - b[i][j]));
    return maxDiff;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv, "params_submanager.input");

    using RotationTypeTag = Properties::TTag::KLPlateTestRotation;
    using DeformationTypeTag = Properties::TTag::KLPlateTestDeformation;
    using CommonTypeTag = Properties::TTag::KLPlateTestCommon;

    using Grid = GetPropType<CommonTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using RotationGridGeometry = GetPropType<RotationTypeTag, Properties::GridGeometry>;
    using DeformationGridGeometry = GetPropType<DeformationTypeTag, Properties::GridGeometry>;
    auto rotationGridGeometry = std::make_shared<RotationGridGeometry>(leafGridView);
    auto deformationGridGeometry = std::make_shared<DeformationGridGeometry>(leafGridView);

    using CouplingManager = GetPropType<RotationTypeTag, Properties::CouplingManager>;
    using RotationProblem = GetPropType<RotationTypeTag, Properties::Problem>;
    using DeformationProblem = GetPropType<DeformationTypeTag, Properties::Problem>;
    using RotationGridVariables = GetPropType<RotationTypeTag, Properties::GridVariables>;
    using DeformationGridVariables = GetPropType<DeformationTypeTag, Properties::GridVariables>;

    using Traits = MultiDomainTraits<RotationTypeTag, DeformationTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;

    constexpr auto rotationIdx = CouplingManager::rotationIdx;
    constexpr auto deformationIdx = CouplingManager::deformationIdx;

    // a deterministic, spatially-varying solution: the plate PDE is linear so the
    // Jacobian does not depend on it, but the coupling terms evaluated through the
    // coupling manager (rotation, deformationAndPotentials) do depend on its values
    auto makeTestSolution = [&](RotationProblem& rotationProblem, DeformationProblem& deformationProblem)
    {
        SolutionVector x;
        rotationProblem.applyInitialSolution(x[rotationIdx]);
        deformationProblem.applyInitialSolution(x[deformationIdx]);

        for (std::size_t i = 0; i < x[rotationIdx].size(); ++i)
            for (int c = 0; c < x[rotationIdx][i].size(); ++c)
                x[rotationIdx][i][c] = std::sin(0.37*i + 0.61*c + 0.2);

        for (std::size_t i = 0; i < x[deformationIdx].size(); ++i)
            for (int c = 0; c < x[deformationIdx][i].size(); ++c)
                x[deformationIdx][i][c] = std::cos(0.29*i - 0.13*c + 0.5)*1e-3;

        return x;
    };

    ////////////////////////////////////////////////////////////////////
    // path A: the regular coupling manager constructed from the grid geometries
    ////////////////////////////////////////////////////////////////////
    auto couplingManagerA = std::make_shared<CouplingManager>(rotationGridGeometry, deformationGridGeometry);
    auto rotationProblemA = std::make_shared<RotationProblem>(rotationGridGeometry, couplingManagerA);
    auto deformationProblemA = std::make_shared<DeformationProblem>(deformationGridGeometry, couplingManagerA);

    auto xA = makeTestSolution(*rotationProblemA, *deformationProblemA);
    couplingManagerA->init(rotationProblemA, deformationProblemA, xA);

    auto rotationGridVariablesA = std::make_shared<RotationGridVariables>(rotationProblemA, rotationGridGeometry);
    auto deformationGridVariablesA = std::make_shared<DeformationGridVariables>(deformationProblemA, deformationGridGeometry);
    rotationGridVariablesA->init(xA[rotationIdx]);
    deformationGridVariablesA->init(xA[deformationIdx]);

    Assembler assemblerA(
        std::make_tuple(rotationProblemA, deformationProblemA),
        std::make_tuple(rotationGridGeometry, deformationGridGeometry),
        std::make_tuple(rotationGridVariablesA, deformationGridVariablesA),
        couplingManagerA
    );

    ////////////////////////////////////////////////////////////////////
    // path B: a default-constructed manager used as a sub-manager of a composed
    // multi-domain coupling manager, attached to the solution storage the composed
    // manager owns (a copy of the solution, like the one its updateSolution fills)
    ////////////////////////////////////////////////////////////////////
    auto couplingManagerB = std::make_shared<CouplingManager>();
    couplingManagerB->computeStencils(rotationGridGeometry, deformationGridGeometry);
    auto rotationProblemB = std::make_shared<RotationProblem>(rotationGridGeometry, couplingManagerB);
    auto deformationProblemB = std::make_shared<DeformationProblem>(deformationGridGeometry, couplingManagerB);

    auto xB = makeTestSolution(*rotationProblemB, *deformationProblemB);
    auto composedStorageB = xB;
    typename CouplingManager::SolutionVectorStorage curSolStorageB = std::make_tuple(
        Dune::stackobject_to_shared_ptr(composedStorageB[rotationIdx]),
        Dune::stackobject_to_shared_ptr(composedStorageB[deformationIdx])
    );
    couplingManagerB->init(rotationProblemB, deformationProblemB, curSolStorageB);

    auto rotationGridVariablesB = std::make_shared<RotationGridVariables>(rotationProblemB, rotationGridGeometry);
    auto deformationGridVariablesB = std::make_shared<DeformationGridVariables>(deformationProblemB, deformationGridGeometry);
    rotationGridVariablesB->init(xB[rotationIdx]);
    deformationGridVariablesB->init(xB[deformationIdx]);

    Assembler assemblerB(
        std::make_tuple(rotationProblemB, deformationProblemB),
        std::make_tuple(rotationGridGeometry, deformationGridGeometry),
        std::make_tuple(rotationGridVariablesB, deformationGridVariablesB),
        couplingManagerB
    );

    ////////////////////////////////////////////////////////////////////
    // compare
    ////////////////////////////////////////////////////////////////////
    const double tol = 1e-10;
    bool failed = false;

    auto check = [&](const std::string& name, double diff)
    {
        std::cout << name << ": max abs difference = " << diff << std::endl;
        if (!(diff <= tol))
            failed = true;
    };

    assemblerA.assembleResidual(xA);
    assemblerB.assembleResidual(xB);
    check("Residual rotation", maxAbsDifferenceVector(assemblerA.residual()[rotationIdx], assemblerB.residual()[rotationIdx]));
    check("Residual deformation", maxAbsDifferenceVector(assemblerA.residual()[deformationIdx], assemblerB.residual()[deformationIdx]));

    assemblerA.assembleJacobianAndResidual(xA);
    assemblerB.assembleJacobianAndResidual(xB);
    check("Jacobian rotation-rotation", maxAbsDifferenceMatrix(
        assemblerA.jacobian()[rotationIdx][rotationIdx], assemblerB.jacobian()[rotationIdx][rotationIdx]));
    check("Jacobian rotation-deformation", maxAbsDifferenceMatrix(
        assemblerA.jacobian()[rotationIdx][deformationIdx], assemblerB.jacobian()[rotationIdx][deformationIdx]));
    check("Jacobian deformation-rotation", maxAbsDifferenceMatrix(
        assemblerA.jacobian()[deformationIdx][rotationIdx], assemblerB.jacobian()[deformationIdx][rotationIdx]));
    check("Jacobian deformation-deformation", maxAbsDifferenceMatrix(
        assemblerA.jacobian()[deformationIdx][deformationIdx], assemblerB.jacobian()[deformationIdx][deformationIdx]));

    // the numeric differentiation deflects the coupling context through the attached
    // storage, which has to hold the solution again afterwards
    check("Attached rotation storage", maxAbsDifferenceVector(composedStorageB[rotationIdx], xB[rotationIdx]));
    check("Attached deformation storage", maxAbsDifferenceVector(composedStorageB[deformationIdx], xB[deformationIdx]));

    if (failed)
    {
        std::cerr << "Sub-manager coupling manager path does not match the regular path!" << std::endl;
        return 1;
    }

    std::cout << "Sub-manager coupling manager path matches the regular path." << std::endl;
    return 0;
}
