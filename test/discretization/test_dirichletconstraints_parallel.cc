// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test of makeConstraintsConsistent: process 0 constrains equation 1 (by primary variable 0) at the dofs of its
 *        interior elements and, in separate entries, equation 0 (by primary variable 1) at those left of x = 0.5; the
 *        last process constrains equation 0 at the dofs of its interior elements. Afterwards every process sharing a
 *        dof constrains each equation that some process constrains, exactly once and with the same value, for PQ2 and
 *        PQ3 finite elements (several dofs per edge and element).
 */
#include <config.h>

#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/constraintinfo.hh>
#include <dumux/common/initialize.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/pq2/fegriddiscretization.hh>
#include <dumux/discretization/pq3/fegriddiscretization.hh>
#include <dumux/linear/parallelhelpers.hh>

namespace {

using namespace Dumux;
constexpr int numEq = 2;
using Values = Dune::FieldVector<double, numEq>;
using Constraint = DirichletConstraintData<DirichletConstraintInfo<numEq>, Values, std::size_t>;

template<class Position>
double valueEq1(const Position& x)
{ return 1.0 + x[0] + 2.0*x[1]; }

template<class Position>
double valueEq0(const Position& x)
{ return 3.0 - x[0] + x[1]; }

template<class GridGeometry>
bool check(const GridGeometry& gridGeometry, const char* label)
{
    const auto& gridView = gridGeometry.gridView();
    constexpr int dim = GridGeometry::GridView::dimension;
    const auto numDofs = gridGeometry.numDofs();
    const int rank = gridView.comm().rank();
    const int lastRank = gridView.comm().size() - 1;

    std::vector<Dune::FieldVector<double, dim>> positions(numDofs);
    std::vector<Constraint> constraints;
    // per dof whether this process constrains equation 0 and equation 1
    Dune::BlockVector<Dune::FieldVector<double, numEq>> isConstrainedSomewhere(numDofs);
    isConstrainedSomewhere = 0.0;
    const auto constrain = [&](std::size_t dofIdx, int eqIdx)
    {
        if (isConstrainedSomewhere[dofIdx][eqIdx] > 0.5)
            return;
        isConstrainedSomewhere[dofIdx][eqIdx] = 1.0;
        const int pvIdx = 1 - eqIdx;
        DirichletConstraintInfo<numEq> info;
        info.set(pvIdx, eqIdx);
        Values values(0.0);
        values[pvIdx] = eqIdx == 1 ? valueEq1(positions[dofIdx]) : valueEq0(positions[dofIdx]);
        constraints.push_back(Constraint{std::move(info), std::move(values), dofIdx});
    };

    auto elemDisc = localView(gridGeometry);
    for (const auto& element : elements(gridView))
    {
        elemDisc.bind(element);
        for (const auto& localDof : localDofs(elemDisc))
            positions[localDof.dofIndex()] = ipData(elemDisc, localDof).global();
    }
    for (const auto& element : elements(gridView, Dune::Partitions::interior))
    {
        elemDisc.bind(element);
        for (const auto& localDof : localDofs(elemDisc))
        {
            if (rank == 0)
                constrain(localDof.dofIndex(), 1);
            if ((rank == 0 && element.geometry().center()[0] < 0.5) || (rank == lastRank && rank != 0))
                constrain(localDof.dofIndex(), 0);
        }
    }
    const auto numLocal = constraints.size();

    makeConstraintsConsistent(gridGeometry, constraints);

    // independent reference: the flags summed over the processes sharing a dof
    using DofMapper = std::decay_t<decltype(gridGeometry.dofMapper())>;
    MultiCodimParallelVectorHelper<typename GridGeometry::GridView, DofMapper> vectorHelper(gridView, gridGeometry.dofMapper());
    vectorHelper.makeNonOverlappingConsistent(isConstrainedSomewhere, std::bitset<dim+1>{(1UL << (dim+1)) - 1UL});

    std::vector<std::array<int, numEq>> count(numDofs, std::array<int, numEq>{});
    bool passed = true;
    for (const auto& constraint : constraints)
    {
        const auto dofIdx = constraint.dofIndex();
        const auto& info = constraint.constraintInfo();
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (!info.isConstraintEquation(eqIdx))
                continue;
            ++count[dofIdx][eqIdx];
            const int pvIdx = info.eqToPriVarIndex(eqIdx);
            const double expected = eqIdx == 1 ? valueEq1(positions[dofIdx]) : valueEq0(positions[dofIdx]);
            passed = passed && pvIdx == 1 - eqIdx && std::abs(constraint.values()[pvIdx] - expected) < 1e-12;
        }
    }
    for (std::size_t i = 0; i < numDofs; ++i)
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            passed = passed && count[i][eqIdx] == (isConstrainedSomewhere[i][eqIdx] > 0.5 ? 1 : 0);

    const auto numAdded = gridView.comm().sum(constraints.size() - numLocal);
    passed = gridView.comm().min(static_cast<int>(passed)) == 1;
    if (rank == 0)
        std::cout << label << ": " << numAdded << " constraints added on the other processes, "
                  << (passed ? "consistent" : "NOT consistent") << std::endl;

    // with several processes, the dofs shared with process 0 are constrained on the other processes too
    return passed && (gridView.comm().size() == 1 || numAdded > 0);
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    const Dune::FieldVector<double, 2> upperRight(1.0);
    const Grid grid(upperRight, {8, 8}, std::bitset<2>{}, /*overlap=*/0);
    const auto gridView = grid.leafGridView();

    bool passed = check(Experimental::PQ2FEGridDiscretization<double, Grid::LeafGridView>(gridView), "PQ2 FE");
    passed = check(Experimental::PQ3FEGridDiscretization<double, Grid::LeafGridView>(gridView), "PQ3 FE") && passed;
    return passed ? 0 : 1;
}
